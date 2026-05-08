# AI-Based Molecule Generation Targeting A2A Adenosine Receptor

Comparing two structure-based generative AI tools — **Pocket2Mol** and **DiffSBDD** — 
for de novo drug candidate generation targeting the A2A adenosine receptor (PDB: 4EIY).

This project was completed as part of the *Artificial Intelligence in Drug Discovery* 
course (Spring 2026) at FHNW University of Applied Sciences and Arts Northwestern Switzerland.

---

## Project Overview

The A2A adenosine receptor is a clinically relevant GPCR target involved in 
Parkinson's disease, immuno-oncology, and cardiac imaging. This project applies 
two generative AI architectures to design novel drug-like molecules that fit 
the receptor's binding pocket:

- **Pocket2Mol** — autoregressive, atom-by-atom generation (SE(3)-equivariant GNN)
- **DiffSBDD** — diffusion-based, all-at-once generation (equivariant score network)

Generated molecules were evaluated using Lipinski's Rule of Five, QED score, 
SA score, and Tanimoto similarity against 10 known A2A actives.

---

---

## Known Limitations & Notes

- DiffSBDD output contained fragmented SMILES (inpainting artifact) — 
  largest fragment selected as primary molecule for analysis
- Pocket2Mol run on CPU (Apple M3) — ~11.5 hours for 113 molecules
- DiffSBDD run required 3 separate runs, outputs combined via `combine_sdfs.py`
- 10 reference A2A actives used for Tanimoto (vs. 14 recommended in assignment)

## My Contributions

| Task | Description |
|---|---|
| PDB Preparation | Downloaded 4EIY, removed ligand/waters/fusion protein in PyMOL, calculated pocket center |
| Pocket2Mol Setup | Cloned repo, adapted environment for Mac ARM64 (CPU-only), ran generation |
| Molecule Generation | Generated 113 valid molecules using Pocket2Mol |
| RDKit Analysis | Calculated MW, LogP, HBA, HBD, TPSA, QED, SA score, Tanimoto similarity |
| Comparative Analysis | Built comparison notebook (04) with overlaid histograms, boxplots, summary table |
| Jupyter Notebooks | Authored 01, 02a, 03a, 04 notebooks with annotated code |

---

## Results Summary

| Metric | Pocket2Mol | DiffSBDD |
|---|---|---|
| Molecules generated | 113 | 116 |
| Lipinski compliant | 113 (100%) | 115 (99.1%) |
| Mean QED | 0.721 | 0.602 |
| Mean SA Score | 3.522 | 4.594 |
| Mean Tanimoto | 0.160 | 0.113 |
| Identical to known active | 1 (AZD4635) | 0 |

Both tools generated novel drug-like molecules (mean Tanimoto < 0.3). 
Pocket2Mol produced higher drug-likeness and better synthesizability. 
DiffSBDD explored a broader, more novel chemical space.

---

## Repository Structure

├── 01_pdb_preparation.ipynb          ← PDB download, PyMOL cleaning, pocket center
├── 02a_pocket2mol_setup.ipynb        ← Pocket2Mol installation and generation
├── 02b_diffsbdd_setup.ipynb          ← DiffSBDD installation and generation
├── 03a_pocket2mol_analysis.ipynb     ← RDKit analysis of Pocket2Mol molecules
├── 03b_diffsbdd_molecule_analysis.ipynb  ← RDKit analysis of DiffSBDD molecules
├── 04_pocket2mol_vs_diffsbdd_comparison.ipynb  ← Comparative analysis
├── data/
│   ├── 4eiy.pdb                      ← original PDB structure
│   └── 4eiy_clean.pdb                ← cleaned structure (ligand/waters removed)
├── outputs/
│   ├── pocket2mol/                   ← Pocket2Mol generated molecules
│   └── diffsbdd/                     ← DiffSBDD generated molecules
├── results/
│   ├── pocket2mol_generated.csv      ← calculated properties (113 molecules)
│   └── diffsbdd_generated.csv        ← calculated properties (116 molecules)
└── figures/                          ← all generated plots and visualizations

---

---

## Generated Figures

| Figure | Description |
|---|---|
| `figures/pocket2mol_property_distributions.png` | MW, LogP, QED, HBA, HBD, TPSA histograms |
| `figures/pocket2mol_tanimoto.png` | Tanimoto distribution + nearest active |
| `figures/pocket2mol_example_molecules.png` | 12 example generated molecules |
| `figures/diffsbdd_property_distributions.png` | Property distributions |
| `figures/diffsbdd_tanimoto.png` | Tanimoto distribution + nearest active |
| `figures/diffsbdd_example_molecules.png` | 12 example generated molecules |
| `figures/comparison_property_distributions.png` | Side-by-side overlaid histograms |
| `figures/comparison_lipinski.png` | Lipinski compliance comparison |
| `figures/comparison_tanimoto_boxplot.png` | Novelty boxplot comparison |
| `figures/comparison_example_molecules.png` | Side-by-side molecule examples |


## Methods

### Target Structure
- PDB: **4EIY** (A2A receptor, 1.8 Å resolution, ZM241385-bound, inactive state)
- Binding pocket center: **(-0.42, 8.53, 17.13)**
- Preparation: removed HOH, ZMA, NA, OLA, OLC, Chain B (BRIL fusion protein)

### Molecular Property Analysis
All properties calculated independently using RDKit:

| Property | Tool | Ideal Range |
|---|---|---|
| MW | `Descriptors.MolWt()` | 150–500 Da |
| LogP | `Descriptors.MolLogP()` | −0.4 – 5.6 |
| HBA/HBD | `NumHAcceptors/Donors()` | ≤10 / ≤5 |
| TPSA | `Descriptors.TPSA()` | ≤140 Å² |
| QED | `QED.qed()` | >0.5 |
| SA Score | `sascorer.calculateScore()` | ≤4 |
| Tanimoto | Morgan fingerprints (r=2, 2048 bits) | <0.3 = novel |

---

## Environment Setup

### Analysis environment (RDKit, pandas, matplotlib)
```bash
conda env create -f requirements_analysis_env.txt
conda activate ai_dd_analysis
```

### Pocket2Mol environment
```bash
conda env create -f env_cuda113_APPLE.yml
conda activate Pocket2Mol
```

### DiffSBDD environment
```bash
conda env create -f DiffSBDD-main/environment.yaml
conda activate diffsbdd
```

---

## References

1. Peng et al. Pocket2Mol. ICML 2022. https://arxiv.org/abs/2205.07249
2. Schneuing et al. DiffSBDD. arXiv:2210.13695, 2022.
3. Jaakola et al. A2A crystal structure. Science 322, 1211 (2008).
4. Bickerton et al. QED. Nature Chemistry 4, 90 (2012).
5. Rogers & Hahn. ECFP fingerprints. JCIM 50, 742 (2010).

---

## Author

**Meryem Topcu**  
MSc Data Science in Life Sciences, FHNW  
meryem.topcu@students.fhnw.ch  
[GitHub](https://github.com/meryemtopcu/ai-a2a-molecule-generation)

*Course: Artificial Intelligence in Drug Discovery, Spring 2026*  
*Supervisor: Dr. Davide Sabbadin | Lecturer: Prof. Dr. Enkelejda Miho*