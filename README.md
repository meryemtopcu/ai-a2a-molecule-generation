# AI in Drug Discovery - Group 2

## Project Summary
This repository compares two structure-based molecular generation pipelines on adenosine A2A receptor structure 4EIY:

1. Pocket2Mol sampling around a cleaned receptor pocket.
2. DiffSBDD inpainting-based generation.
3. Post-processing and property comparison, including Tanimoto similarity against known A2A actives.

## Deliverables
### Task 1: Problem Definition
- Output: 2-page report (PDF/Word, APA references).
- Topics: A2A receptor background, motivation for new ligands, method overview, alternatives.

### Task 2: Method Execution
- Output: annotated Python workflow and generated molecule sets.
- Steps: prepare 4EIY, run Pocket2Mol (>=100 molecules), run DiffSBDD (>=100 molecules), compute properties, compare against known A2A actives.

### Task 3: Results Presentation
- Output: 20 min talk + 5 min Q&A.
- Content: intro, methods, results, comparison, conclusions.

## Repository Map
- Main property script: [compute_properties.py](compute_properties.py)
- Pocket2Mol entrypoint: [Pocket2Mol/sample_for_pdb.py](Pocket2Mol/sample_for_pdb.py)
- Pocket2Mol config: [Pocket2Mol/configs/sample_for_pdb.yml](Pocket2Mol/configs/sample_for_pdb.yml)
- Pocket2Mol checkpoint location: [Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt](Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt)
- DiffSBDD entrypoint: [DiffSBDD-main/inpaint.py](DiffSBDD-main/inpaint.py)
- DiffSBDD environment: [DiffSBDD-main/environment.yaml](DiffSBDD-main/environment.yaml)
- DiffSBDD checkpoint location: [DiffSBDD-main/checkpoints/crossdocked_fullatom_cond.ckpt](DiffSBDD-main/checkpoints/crossdocked_fullatom_cond.ckpt)
- A2A actives SDF: [outputs2/a2a_drugs/ADORA2A-world.sdf](outputs2/a2a_drugs/ADORA2A-world.sdf)
- 2D visualization notebook: [visualise_2D_sdf.ipynb](visualise_2D_sdf.ipynb)

## External Resources
- PDB entry 4EIY: https://www.rcsb.org/structure/4EIY
- Pocket2Mol paper: https://arxiv.org/abs/2205.07249
- DiffSBDD paper: https://arxiv.org/abs/2210.13695
- DiffSBDD pretrained checkpoints: https://zenodo.org/record/8183747
- A2A actives subset (ZINC): https://zinc15.docking.org/genes/ADORA2A/substances/subsets/world/
- A2A world SDF (all molecules): https://zinc15.docking.org/genes/ADORA2A/substances/subsets/world.sdf?count=all
- PyMOL: https://pymol.org/
- RDKit: https://www.rdkit.org/

## Environment Setup
### Pocket2Mol environment
Use the Apple-compatible environment file:
- [env_cuda113_APPLE.yml](env_cuda113_APPLE.yml)

```bash
conda env create -f env_cuda113_APPLE.yml
conda activate Pocket2Mol
```

### DiffSBDD environment
```bash
conda env create -f DiffSBDD-main/environment.yaml
conda activate diffsbdd
```

## 4EIY Preparation in PyMOL
Use the original structure, then create cleaned protein coordinates.

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

Pocket center used in this project:
- x = -0.42
- y = 8.53
- z = 17.13

Cleaned structure used by scripts:
- [Pocket2Mol/data/4eiy_clean.pdb](Pocket2Mol/data/4eiy_clean.pdb)

## Run Pocket2Mol
Run from repository root after placing checkpoint in [Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt](Pocket2Mol/ckpt/pretrained_Pocket2Mol.pt).

```bash
python Pocket2Mol/sample_for_pdb.py \
  --pdb_path Pocket2Mol/data/4eiy_clean.pdb \
  --config Pocket2Mol/configs/sample_for_pdb.yml \
  --center " -0.42,8.53,17.13"
```

Notes:
- Leading space before negative x is intentional.
- To customize outputs, pass `--outdir /path/to/folder`.

## Run DiffSBDD
```bash
conda run -n diffsbdd python DiffSBDD-main/inpaint.py \
  DiffSBDD-main/checkpoints/crossdocked_fullatom_cond.ckpt \
  --pdbfile Pocket2Mol/data/4eiy_clean.pdb \
  --outfile outputs2/4eiy_inpaint.sdf \
  --ref_ligand DiffSBDD-main/example/5ndu_C_8V2.sdf \
  --fix_atoms DiffSBDD-main/example/fragments.sdf \
  --center pocket \
  --n_samples 100
```

## Compute and Compare Properties
Script:
- [compute_properties.py](compute_properties.py)

Current workflow in script:
1. Loads Pocket2Mol outputs from a PATH1 SDF folder.
2. Loads DiffSBDD outputs from a PATH2 multi-molecule SDF.
3. Loads known A2A actives from [outputs2/a2a_drugs/ADORA2A-world.sdf](outputs2/a2a_drugs/ADORA2A-world.sdf).
4. Validates molecules and reports invalid counts per path.
5. Computes:
- MolecularWeight
- LogP
- NumHAcceptors
- NumHDonors
- TPSA
- QED
- LipinskiRuleOfFive
- TanimotoSimilarityA2A

Run:
```bash
python compute_properties.py
```

Output JSON default:
- [outputs2/compute_properties.json](outputs2/compute_properties.json)

## Visualize Generated Molecules
- 2D grid/inspection notebook: [visualise_2D_sdf.ipynb](visualise_2D_sdf.ipynb)
- Optional PyMOL 3D inspection of receptor + ligands.

## Notes on macOS Adaptation
- Pocket2Mol code in this repository includes CPU-friendly changes for non-CUDA execution.
- DiffSBDD environment has been adapted for macOS-compatible dependency resolution.

## Recommended Next Improvements
1. Move hardcoded paths in [compute_properties.py](compute_properties.py) to CLI arguments.
2. Add histogram/heatmap generation script for report-ready figures.
3. Add mean/top-k A2A similarity summary per method in output JSON.