import pandas as pd
from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers, inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from src.curator_2d import _write_sdf

try:
    from dimorphite_dl import DimorphiteDL
    HAS_DIMORPHITE = True
except ImportError:
    HAS_DIMORPHITE = False


def step7_physiological_ionization(input_sdf: str, output_sdf: str, output_csv: str, ph_min: float = 7.2, ph_max: float = 7.6) -> pd.DataFrame:
    """Step 7: Assigns biological microspecies at pH 7.4 ± 0.2 via Dimorphite-DL (with Uncharger fallback)."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    uncharger = rdMolStandardize.Uncharger()
    dimorphite = DimorphiteDL(min_ph=ph_min, max_ph=ph_max, pka_precision=1.0) if HAS_DIMORPHITE else None

    ionized_mols, records = [], []
    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        ik = mol.GetProp("InChIKey_Hash") if mol.HasProp("InChIKey_Hash") else "N/A"
        base_smi = Chem.MolToSmiles(mol, isomericSmiles=True)

        if dimorphite:
            protonated_smiles = dimorphite.protonate(base_smi)
        else:
            protonated_smiles = [Chem.MolToSmiles(uncharger.uncharge(mol), isomericSmiles=True)]

        for p_idx, p_smi in enumerate(protonated_smiles):
            p_mol = Chem.MolFromSmiles(p_smi)
            if p_mol is None:
                continue
            micro_name = f"{name}_pH_{p_idx+1}" if len(protonated_smiles) > 1 else name
            p_mol.SetProp("_Name", micro_name)
            p_mol.SetProp("InChIKey_Hash", inchi.MolToInchiKey(p_mol))
            ionized_mols.append(p_mol)
            records.append({
                "Compound_ID": micro_name,
                "InChIKey_Hash": p_mol.GetProp("InChIKey_Hash"),
                "Ionized_SMILES": p_smi,
                "Formal_Charge": Chem.GetFormalCharge(p_mol)
            })

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(ionized_mols, output_sdf)
    return df


def step8_enumerate_stereoisomers(input_sdf: str, output_sdf: str, output_csv: str, max_isomers: int = 32) -> pd.DataFrame:
    """Step 8: Combinatorially enumerates unassigned chiral centers and E/Z double bonds."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    enum_mols, records = [], []

    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        generator = EnumerateStereoisomers.EnumerateStereoisomers(mol)

        isomers = []
        for iso in generator:
            isomers.append(iso)
            if len(isomers) >= max_isomers:
                break

        for iso_idx, iso_mol in enumerate(isomers):
            iso_name = f"{name}_stereo_{iso_idx+1}"
            iso_ik = inchi.MolToInchiKey(iso_mol)
            iso_mol.SetProp("_Name", iso_name)
            iso_mol.SetProp("InChIKey_Hash", iso_ik)
            enum_mols.append(iso_mol)
            records.append({
                "Parent_ID": name,
                "Isomer_ID": iso_name,
                "Isomer_InChIKey": iso_ik,
                "Isomer_SMILES": Chem.MolToSmiles(iso_mol, isomericSmiles=True)
            })

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(enum_mols, output_sdf)
    return df
