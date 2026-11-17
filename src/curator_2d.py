import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, FilterCatalog, inchi
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalogParams


def step1_deduplicate(input_sdf: str, output_sdf: str, output_csv: str) -> pd.DataFrame:
    """Step 1: Canonicalizes structures and deduplicates via deterministic InChIKey hashes."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    seen_keys = set()
    unique_mols, records = [], []

    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") and mol.GetProp("_Name").strip() else f"Ligand_{idx+1}"
        try:
            ik = inchi.MolToInchiKey(mol)
            if not ik:
                raise ValueError("Empty InChIKey")
        except Exception:
            ik = inchi.MolToInchiKey(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))

        if ik in seen_keys:
            continue
        seen_keys.add(ik)

        mol.SetProp("_Name", name)
        mol.SetProp("InChIKey_Hash", ik)
        unique_mols.append(mol)
        records.append({"Compound_ID": name, "InChIKey_Hash": ik, "SMILES": Chem.MolToSmiles(mol, isomericSmiles=True)})

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(unique_mols, output_sdf)
    return df


def step2_sanitize(input_sdf: str, output_sdf: str, output_csv: str) -> pd.DataFrame:
    """Step 2: Audits chemical valency and standardizes aromaticity models."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    valid_mols, records = [], []

    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        ik = mol.GetProp("InChIKey_Hash") if mol.HasProp("InChIKey_Hash") else "N/A"
        try:
            san_mol = Chem.Mol(mol)
            Chem.SanitizeMol(san_mol)
            san_mol.SetProp("_Name", name)
            san_mol.SetProp("InChIKey_Hash", ik)
            valid_mols.append(san_mol)
            records.append({"Compound_ID": name, "InChIKey_Hash": ik, "Sanitized_SMILES": Chem.MolToSmiles(san_mol, isomericSmiles=True)})
        except Exception:
            continue

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(valid_mols, output_sdf)
    return df


def step3_desalt(input_sdf: str, output_sdf: str, output_csv: str) -> pd.DataFrame:
    """Step 3: Strips counter-ions and residual crystallization solvents."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    remover = SaltRemover()
    clean_mols, records = [], []

    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        ik = mol.GetProp("InChIKey_Hash") if mol.HasProp("InChIKey_Hash") else "N/A"
        try:
            desalted = remover.StripMol(mol, dontRemoveEverything=True)
            stripped = desalted.GetNumAtoms() < mol.GetNumAtoms()
            desalted.SetProp("_Name", name)
            desalted.SetProp("InChIKey_Hash", ik)
            clean_mols.append(desalted)
            records.append({"Compound_ID": name, "InChIKey_Hash": ik, "Desalted_SMILES": Chem.MolToSmiles(desalted, isomericSmiles=True), "Salt_Stripped": stripped})
        except Exception:
            continue

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(clean_mols, output_sdf)
    return df


def step4_standardize(input_sdf: str, output_sdf: str, output_csv: str) -> pd.DataFrame:
    """Step 4: Normalizes functional groups, mesomeric forms, and tautomeric states."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    tautomer_canon = rdMolStandardize.TautomerEnumerator()
    std_mols, records = [], []

    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        ik = mol.GetProp("InChIKey_Hash") if mol.HasProp("InChIKey_Hash") else "N/A"
        try:
            clean_mol = rdMolStandardize.Cleanup(mol)
            clean_mol = tautomer_canon.Canonicalize(clean_mol)
            clean_mol.SetProp("_Name", name)
            clean_mol.SetProp("InChIKey_Hash", ik)
            std_mols.append(clean_mol)
            records.append({"Compound_ID": name, "InChIKey_Hash": ik, "Standardized_SMILES": Chem.MolToSmiles(clean_mol, isomericSmiles=True)})
        except Exception:
            continue

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(std_mols, output_sdf)
    return df


def step5_bioavailability_filter(input_sdf: str, output_sdf: str, output_csv: str, mode: str = "bro5") -> pd.DataFrame:
    """Step 5: Filters chemical space using strict Ro5/Veber or marine-tolerant bRo5 cutoffs."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    passed_mols, records = [], []

    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        ik = mol.GetProp("InChIKey_Hash") if mol.HasProp("InChIKey_Hash") else "N/A"

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotb = Lipinski.NumRotatableBonds(mol)

        if mode == "ro5":
            violations = sum([mw > 500, logp > 5.0, hbd > 5, hba > 10])
            passed = (violations <= 1) and (rotb <= 10) and (tpsa <= 140.0)
        else:  # Beyond-Rule-of-5 (bRo5) marine natural product regime
            passed = (mw <= 800.0) and (-2.0 <= logp <= 6.5) and (hbd <= 8) and (hba <= 12) and (tpsa <= 180.0)

        if passed:
            for prop, val in [("MW", f"{mw:.2f}"), ("LogP", f"{logp:.2f}"), ("TPSA", f"{tpsa:.2f}"), ("RotB", str(rotb))]:
                mol.SetProp(prop, val)
            passed_mols.append(mol)
            records.append({"Compound_ID": name, "InChIKey_Hash": ik, "MW": round(mw, 2), "LogP": round(logp, 2), "HBD": hbd, "HBA": hba, "TPSA": round(tpsa, 2), "RotB": rotb})

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(passed_mols, output_sdf)
    return df


def step6_pains_filter(input_sdf: str, output_sdf: str, output_csv: str) -> pd.DataFrame:
    """Step 6: Purges Pan-Assay Interference Compounds (PAINS A, B, C)."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)

    clean_mols, records = [], []
    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{idx+1}"
        ik = mol.GetProp("InChIKey_Hash") if mol.HasProp("InChIKey_Hash") else "N/A"

        if not catalog.HasMatch(mol):
            clean_mols.append(mol)
            records.append({"Compound_ID": name, "InChIKey_Hash": ik, "PAINS_Status": "PASS"})

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    _write_sdf(clean_mols, output_sdf)
    return df


def _write_sdf(mols, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    writer = Chem.SDWriter(path)
    for m in mols:
        writer.write(m)
    writer.close()
