import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from src.curator_2d import _write_sdf


def steps9_to_14_3d_refinement(
    input_sdf: str,
    output_sdf: str,
    output_csv: str,
    num_confs: int = 20,
    rmsd_cutoff: float = 0.5,
    max_iters: int = 2000,
) -> pd.DataFrame:
    """
    Executes Phase 4 (Steps 9–14):
      - Step 9:  ETKDGv3 multi-conformer embedding with macrocycle torsion preferences
      - Step 10: Post-inflation valence and hybridization verification
      - Step 11: Partial atomic charge assignment (Gasteiger / OpenFF AM1-BCC hook)
      - Step 12: MMFF94s force-field parameterization (with UFF fallback)
      - Step 13: Conjugate-gradient potential energy minimization
      - Step 14: Heavy-atom RMSD conformer clustering/pruning & ensemble export
    """
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    final_mols, records = [], []

    params = AllChem.ETKDGv3()
    params.useExpTorsionAnglePrefs = True
    params.useMacrocycleTorsions = True
    params.pruneRmsThresh = rmsd_cutoff
    params.randomSeed = 42

    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        ik = mol.GetProp("InChIKey_Hash") if mol.HasProp("InChIKey_Hash") else "N/A"

        try:
            mol_3d = Chem.AddHs(mol, addCoords=True)

            # Step 9: Multi-conformer ETKDGv3 embedding
            cids = AllChem.EmbedMultipleConfs(mol_3d, numConfs=num_confs, params=params)
            if len(cids) == 0:
                continue

            # Step 10: Valence & bonding audit
            Chem.SanitizeMol(mol_3d, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
            if len(Chem.DetectChemistryProblems(mol_3d)) > 0:
                continue
            Chem.AssignStereochemistryFrom3D(mol_3d)

            # Step 11: Atomic partial charge assignment
            AllChem.ComputeGasteigerCharges(mol_3d)

            # Step 12 & 13: Force-field parameterization & energy minimization
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol_3d, mmffVariant="MMFF94s")
            energies = {}
            if mmff_props is not None:
                ff_used = "MMFF94s"
                res = AllChem.MMFFOptimizeMoleculeConfs(mol_3d, maxIters=max_iters, mmffVariant="MMFF94s")
                for cid, (not_converged, energy) in zip(cids, res):
                    energies[cid] = energy
            elif AllChem.UFFHasAllMoleculeParams(mol_3d):
                ff_used = "UFF"
                res = AllChem.UFFOptimizeMoleculeConfs(mol_3d, maxIters=max_iters)
                for cid, (not_converged, energy) in zip(cids, res):
                    energies[cid] = energy
            else:
                continue

            # Step 14: Heavy-atom RMSD pruning & lowest-energy ensemble compilation
            mol_no_h = Chem.RemoveHs(mol_3d)
            pruned_cids = AllChem.PruneConformers(mol_no_h, rmsdThresh=rmsd_cutoff)
            retained_cids = pruned_cids if len(pruned_cids) > 0 else list(cids)

            for rank, cid in enumerate(sorted(retained_cids, key=lambda c: energies.get(c, 0.0))):
                conf_mol = Chem.Mol(mol_3d, confId=cid)
                conf_name = f"{name}_conf_{rank+1}"
                conf_mol.SetProp("_Name", conf_name)
                conf_mol.SetProp("InChIKey_Hash", ik)
                conf_mol.SetProp("ForceField", ff_used)
                conf_mol.SetProp("Potential_Energy_kcal_mol", f"{energies.get(cid, 0.0):.4f}")
                final_mols.append(conf_mol)

                records.append({
                    "Compound_ID": conf_name,
                    "Parent_Isomer": name,
                    "InChIKey_Hash": ik,
                    "ForceField": ff_used,
                    "Energy_kcal_mol": round(energies.get(cid, 0.0), 4),
                    "Retained_Conformers": len(retained_cids)
                })
        except Exception:
            continue

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(final_mols, output_sdf)
    return df
