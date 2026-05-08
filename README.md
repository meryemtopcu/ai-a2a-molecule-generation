# AI in Drug Discovery - Group 2

## Project Summary
This repository compares two structure-based molecular generation pipelines on adenosine A2A receptor structure 4EIY:

1. Pocket2Mol sampling around a cleaned receptor pocket.
2. DiffSBDD inpainting-based generation.
3. Post-processing and property comparison, including Tanimoto similarity against known A2A actives.
4. Automatic merging of split SDF outputs for unified analysis.

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
- SDF merging utility: [combine_sdfs.py](combine_sdfs.py)
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
Use `compute_properties.py` to validate generated molecules, compute standard descriptors, and record similarity to known A2A actives.

What it does (high level)
- Loads generated molecules (Pocket2Mol one-file-per-molecule and DiffSBDD multi-molecule SDFs).
- Loads an A2A reference set (ZINC world subset at `outputs2/a2a_drugs/ADORA2A-world.sdf`).
- Validates molecules, computes descriptors, and saves per-molecule CSVs.

Default inputs
- Pocket2Mol samples (default PATH1): `outputs2/sample_for_pdb_4EIY_clean.pdb_2026_05_02__11_02_39/SDF`
- DiffSBDD combined SDF (default PATH2): `outputs2/combined_diffsbdd.sdf`
- A2A reference SDF (default `--a2a`): `outputs2/a2a_drugs/ADORA2A-world.sdf` (downloaded from ZINC15)

Outputs
- `results/pocket2mol_generated.csv` (if PATH1 used)
- `results/diffsbdd_generated.csv`

Descriptors computed (per valid molecule)
- `smiles`, `MW`, `LogP`, `HBA`, `HBD`, `TPSA`, `RotBonds`, `QED`, `SA_score`, `Lipinski`
- Similarity: `max_tanimoto` (maximum Tanimoto vs. A2A reference; Morgan radius=2, 2048 bits)

Quick commands
- Default run:
```bash
python compute_properties.py
```
- Custom inputs:
```bash
python compute_properties.py --path1 /path/to/pocket2mol_dir --path2 /path/to/diffsbdd.sdf --a2a /path/to/ADORA2A-world.sdf
```

Notes
- Fingerprints: Morgan (r=2, nBits=2048).
- SA score is optional (requires `sascorer`); if unavailable the field is `None`.
- Use `--plots-only` to regenerate plots from existing CSVs without recomputing descriptors.


## Similarity Analysis vs. Known Actives
Use the direct SDF-based similarity analysis script to compare Pocket2Mol and DiffSBDD against the 14 ChEMBL reference actives in `ChEMBL reference molecules`.

Run the full analysis:
```bash
python similarity_analysis.py
```

This reads:
- Pocket2Mol SDFs from `outputs2/sample_for_pdb_4EIY_clean.pdb_2026_05_02__11_02_39/SDF`
- DiffSBDD SDFs from `outputs2/combined_diffsbdd.sdf`
- 14 ChEMBL reference ligands from `ChEMBL reference molecules/DOWNLOAD-muQF1iIRyjhMqiuLdRuitIuq-3YXqGwME0QTLPrZB4E_eq_.csv`

It writes the combined comparison figure to `results/similarity_vs_known_actives.png`.

Run with Pocket2Mol analysis (provide PATH1):
```bash
python compute_properties.py \
  --path1 /path/to/pocket2mol/samples
```
Analyzes Pocket2Mol output from PATH1. If PATH2 is also provided, both pipelines are processed independently.

Run with fully custom paths:
```bash
python compute_properties.py \
  --path1 /path/to/pocket2mol/samples \
  --path2 /path/to/4eiy_inpaint.sdf \
  --a2a /path/to/ADORA2A-world.sdf
```

CLI arguments (all optional):
- `--path1`: Directory of single-molecule SDF files from Pocket2Mol; if omitted, PATH1 is skipped unless both paths are omitted, in which case the default PATH1 location is used.
- `--path2`: Multi-molecule SDF file from DiffSBDD; if omitted, PATH2 is skipped unless both paths are omitted, in which case the default PATH2 location is used.
- `--a2a`: A2A actives reference SDF (default: `outputs2/a2a_drugs/ADORA2A-world.sdf`)
- `--append`: Append PATH2 results to `diffsbdd_generated.csv` instead of overwriting it.
- `--plots-only`: Recompute similarity plots without rewriting the CSV outputs.

Output CSV files:
- [results/pocket2mol_generated.csv](results/pocket2mol_generated.csv) (if PATH1 provided)
- [results/diffsbdd_generated.csv](results/diffsbdd_generated.csv) (always, unless you run `--plots-only`)

CSV columns: `smiles`, `valid`, `MW`, `LogP`, `HBA`, `HBD`, `TPSA`, `RotBonds`, `QED`, `SA_score`, `Lipinski`, `max_tanimoto`

## `compute_properties.py` Function Reference

Property and validation helpers:
- `_compute_all_properties(mol)`: Computes the full descriptor bundle used everywhere else, including MW, LogP, HBA, HBD, TPSA, rotatable bonds, QED, SA score, and the Lipinski rule count.
- `is_valid_molecule(mol)`: Runs the descriptor bundle inside a safety check and returns whether the molecule is processable.
- `compute_sa_score(mol)`: Calls the SA score helper when available and returns `None` if the scorer cannot be loaded or fails.
- `get_smiles(mol)`: Converts a molecule into canonical SMILES.
- `compute_properties(mol)`: Public wrapper around the internal descriptor bundle for validated molecules.

Similarity and fingerprint helpers:
- `compute_tanimoto_similarity(mol, reference_data)`: Computes the maximum Tanimoto similarity between one generated molecule and the A2A reference set.
- `get_morgan_fingerprint(mol)`: Sanitizes a copy of the molecule and builds the Morgan fingerprint used for similarity calculations.
- `load_a2a_reference_fingerprints(a2a_sdf_file)`: Loads the A2A reference SDF and returns fingerprints plus reference names.

Input loading helpers:
- `load_molecules_from_path1(path1)`: Loads one-molecule-per-file Pocket2Mol SDFs from a directory.
- `load_molecules_from_path2(path2)`: Loads all molecules from a single multi-molecule DiffSBDD SDF.
- `filter_duplicates(molecules)`: Removes duplicate molecules based on canonical SMILES.

CSV and comparison helpers:
- `save_properties_to_csv(molecules_data, output_file, append=False)`: Writes the property table to CSV and supports append mode for PATH2.
- `compare_properties(properties1, properties2)`: Compares two property dictionaries and returns only the differing entries.
- `compare_property_lists(properties1, properties2)`: Compares paired molecule property lists entry by entry.

Similarity reporting and plotting helpers:
- `load_similarity_rows(csv_file)`: Reads `max_tanimoto` values back from a CSV.
- `compute_similarity_matrix(molecules, reference_data)`: Builds a full generated-vs-reference Tanimoto matrix for heatmap plotting.
- `plot_similarity_comparison(...)`: Produces the max-Tanimoto distribution plot, the binned tool-comparison heatmap, and the full tool-vs-reference heatmaps.

Workflow entry point:
- `main()`: Parses CLI arguments, loads inputs, validates molecules, computes properties, writes CSVs, and triggers plot generation.
  - *Note*: If `PATH2` is a directory, it automatically triggers `combine_sdfs` before processing.

## Visualize Generated Molecules
- 2D grid/inspection notebook: [visualise_2D_sdf.ipynb](visualise_2D_sdf.ipynb)
- Optional PyMOL 3D inspection of receptor + ligands.

## Notes on macOS Adaptation
- Pocket2Mol code in this repository includes CPU-friendly changes for non-CUDA execution.
- DiffSBDD environment has been adapted for macOS-compatible dependency resolution.
- `compute_properties.py` is the active analysis entry point; the two standalone RDKit helper drafts above are not imported by the current workflow.

## Recommended Next Improvements