# AI in drug discovery - Group2

# SUMMARY
## Task 1: Problem Definition ✅ (attached docx file)
Output: 2-page report (PDF/Word) with APA references
What to include:
- What is the A2A receptor and why is it important?
- Why do we need new molecules?
- What are Pocket2Mol and DiffSBDD?
- What other methods exist?
  
## Task 2: Method Training 🔧
Output: Annotated source code (Python)
What to do:
- Download PDB 4EIY crystal structure
- Clean the structure (remove ligand, water molecules)
- Run Pocket2Mol → generate ≥100 molecules
- Run DiffSBDD → generate ≥100 molecules
- Compute molecular properties (MW, LogP, QED, etc.)
- Tanimoto similarity vs. known A2A actives
- Create histograms and heatmaps

## Task 3: Results Summary 🎤
Output: Presentation in front of class (20 min + 5 min Q&A)
Slides should cover:
- Introduction (A2A biology)
- Methods (Pocket2Mol vs DiffSBDD)
- Results (generated molecules, property analysis)
- Discussion (which tool is better and why?)
- Conclusion & outlook

## Getting started
Install pymol
'''clean the file with
load 4EIY.pdb
remove solvent
remove inorganic
select prot, polymer.protein
create clean_protein, prot
save 4EIY_clean.pdb, clean_protein'''

### How to find the binding pocket
The Pocket2Mol and DiffSBDD runs use the center of the co-crystallized ligand ZM241385 before it is removed from the structure.

Workflow in PyMOL:

```pymol
load 4EIY.pdb
remove solvent
remove inorganic
remove resn OLA+OLC
remove chain B
select prot, polymer.protein
create clean_protein, prot
save 4EIY_clean.pdb, clean_protein
```

The pocket center used in this project is:

```text
x = -0.42, y = 8.53, z = 17.13
```

These coordinates correspond to the ligand position and can be passed directly to `sample_for_pdb.py`.

## Install pocket2mol
There are several architectures, the conda.env is for OSARM64 and does not work with Apple M Processors. Use the file [env_cuda113_APPLE.yml](env_cuda113_APPLE.yml)

### Install via conda yaml file (cuda 11.3)
'''conda env create -f env_cuda113_APPLE.yml
conda activate Pocket2Mol'''

### macOS-specific Pocket2Mol changes
This repository has a few changes so Pocket2Mol can run on macOS without CUDA:

```text
- checkpoint paths are resolved relative to Pocket2Mol/ so you can run from the repo root
- CUDA is optional; the sampling scripts fall back to CPU when torch.cuda is unavailable
- torch-cluster / torch-scatter calls have pure PyTorch fallbacks for KNN and scatter operations
```

The pretrained checkpoint is expected at [Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt](Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt).

## PDB Structure Preparation — 4EIY

### Overview
The crystal structure of the A2A adenosine receptor (PDB: 4EIY) was 
downloaded from the RCSB Protein Data Bank. The structure was solved 
at 1.8 Å resolution and contains the co-crystallized antagonist ZM241385.

## Cleaning Steps (performed in PyMOL)
The following components were removed to obtain a clean receptor structure:

- **HOH** — water molecules (185 atoms removed)
- **ZMA** — co-crystallized ligand ZM241385 (25 atoms removed)
- **NA** — sodium ion (1 atom removed)
- **OLA, OLC** — lipid molecules from membrane mimetic
- **Chain B** — BRIL fusion protein used for crystallization

## Output
Clean structure saved as: `data/4eiy_clean.pdb`
Original structure retained as: `data/4eiy.pdb`

## Binding Pocket Center
The center of the binding pocket was determined from the ZM241385 
ligand position prior to removal:
- Coordinates: **(x=-0.42, y=8.53, z=17.13)**

### How to run Pocket2Mol
Run the sampling command from the repository root after placing the checkpoint in [Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt](Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt):

```bash
python Pocket2Mol/sample_for_pdb.py \
	--pdb_path Pocket2Mol/data/4EIY_clean.pdb \
	--config Pocket2Mol/configs/sample_for_pdb.yml \
	--center " -0.42,8.53,17.13"
```

The leading space before the first coordinate is required because the x-value is negative.
- These coordinates will be used as input for Pocket2Mol and DiffSBDD


