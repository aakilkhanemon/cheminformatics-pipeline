# MarDockPrep: Cheminformatics Pipeline for Marine Natural Product Libraries

![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Python 3.12+](https://img.shields.io/badge/Python-3.12+-blue.svg)
![RDKit 2026.03](https://img.shields.io/badge/RDKit-2026.03-orange.svg)

![MarDockPrep Workflow](assets/workflow_architecture.png)

Official implementation of the 14-step RDKit-accelerated small-molecule curation, stereochemical enumeration, and 3D force-field refinement pipeline from:  
**"MarDockPrep: An Automated 14-Step Cheminformatics Pipeline for Preparing Complex Marine Chemical Libraries"**  
*Mohammed Emon et al.*

---

## Algorithmic Architecture & 14-Step Decision Logic

Marine secondary metabolites frequently exhibit dense stereocenters, macrocyclic backbones, heavy halogenation, and Beyond-Rule-of-5 (bRo5) physicochemical profiles. **MarDockPrep** automates library preparation across three modular computational stages:

1. **2D Structural Curation & Property Filtering (`src/curator_2d.py`)**
   * **Step 1 (Deduplication):** Deterministic `InChIKey` hashing (`rdkit.Chem.inchi`) to eliminate redundant entries.
   * **Step 2 (Sanitization):** Strict valence and aromaticity verification (`Chem.SanitizeMol`).
   * **Step 3 (Desalting):** Removal of counter-ions and residual solvents (`SaltRemover.StripMol`).
   * **Step 4 (Standardization):** Functional group and tautomer canonicalization (`rdMolStandardize`).
   * **Step 5 (Bioavailability Filtering):** Configurable Lipinski/Veber (`--filter-mode ro5`) or marine-tolerant bRo5 (`--filter-mode bro5`) gating.
   * **Step 6 (PAINS Exclusion):** Substructure filtering against `FilterCatalog` PAINS A/B/C catalogs.

2. **Ionization & Stereochemical Expansion (`src/expander_chem.py`)**
   * **Step 7 (Physiological Protonation):** Microspecies generation at $\text{pH } 7.4 \pm 0.2$ via `Dimorphite-DL`.
   * **Step 8 (Stereochemical Enumeration):** Explicit expansion of undefined chiral centers and $E/Z$ bonds (`EnumerateStereoisomers`, capped via `--max-isomers`).

3. **3D Coordinate Generation & Ensemble Refinement (`src/builder_3d.py`)**
   * **Step 9 (3D Embedding):** Distance geometry embedding with experimental macrocycle torsions (`ETKDGv3`).
   * **Step 10 (Valence & Hybridization Audit):** Post-inflation geometry verification (`DetectChemistryProblems`).
   * **Step 11 (Partial Charges):** Atomic partial charge assignment.
   * **Step 12–13 (Force-Field Parameterization & Minimization):** `MMFF94s` (with `UFF` fallback) conjugate-gradient relaxation.
   * **Step 14 (Ensemble Pruning):** Heavy-atom RMSD clustering (`AllChem.PruneConformers`) and multi-model `.sdf` compilation.

---

## Quickstart & CLI Execution

```bash
# 1. Clone repository and create environment
git clone [https://github.com/aakilkhanemon/mardockprep.git](https://github.com/aakilkhanemon/mardockprep.git)
cd mardockprep
conda env create -f environment.yml
conda activate mardockprep

# 2. Execute end-to-end 14-step preparation pipeline
python -m src.pipeline_runner \
  --input data/raw/FDAAI_20_Meta.sdf \
  --outdir data/processed \
  --filter-mode bro5 \
  --max-isomers 32 \
  --num-confs 20 \
  --rmsd-cutoff 0.5
```

---

## Directory Layout

```text
mardockprep/
├── assets/
│   └── workflow_architecture.png     # Workflow and 14-step architecture diagrams
├── data/
│   ├── raw/                          # Input raw marine SMILES / SDF datasets
│   └── processed/                    # Step-wise intermediate & docking-ready 3D SDF files
├── notebooks/
│   └── mardockprep_colab.ipynb       # Interactive Google Colab notebook
├── src/
│   ├── __init__.py                   # Package initialization
│   ├── curator_2d.py                 # Steps 1–6: 2D filtration, desalting, standardization & PAINS
│   ├── expander_chem.py              # Steps 7–8: Physiological ionization & stereoisomer enumeration
│   ├── builder_3d.py                 # Steps 9–14: ETKDGv3 embedding, charge/forcefield & RMSD pruning
│   └── pipeline_runner.py            # Master CLI orchestrator (`python -m src.pipeline_runner`)
├── results/
│   ├── figures/                      # Exported 600 DPI PNG/TIFF figures
│   └── tables/                       # Step-by-step audit CSVs (Table 1 & Table 2)
├── CITATION.cff                      # Academic citation metadata
├── environment.yml                   # Reproducible Conda environment
├── LICENSE                           # MIT License
└── README.md                         # Repository documentation
```

---

## Citation

If you use **MarDockPrep** in your research, please cite:

```bibtex
@software{emon2026mardockprep,
  author       = {Emon, Mohammed},
  title        = {{MarDockPrep: An Automated 14-Step Cheminformatics Pipeline for Preparing Complex Marine Chemical Libraries}},
  year         = {2026},
  version      = {1.0.0},
  url          = {[https://github.com/aakilkhanemon/mardockprep](https://github.com/aakilkhanemon/mardockprep)}
}
```
