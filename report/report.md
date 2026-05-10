# Written Report Structure

Target length: 10-20 pages

## Abstract
The report compares Pocket2Mol and DiffSBDD, for ligand design against the A2A adenosine receptor (PDB 4EIY). Both methods were run on the same prepared receptor pocket and evaluated with a shared post-processing pipeline to enable a fair comparison. Generated molecules were assessed using RDKit descriptors (MW, LogP, HBA, HBD, TPSA, rotatable bonds), QED, synthetic accessibility (SA score), Lipinski compliance, and Tanimoto similarity to known A2A reference actives using Morgan fingerprints (radius 2, 2048 bits).

Pocket2Mol produced 113 valid molecules and DiffSBDD produced 116 valid molecules. Lipinski compliance was high for both tools (Pocket2Mol: 100%, DiffSBDD: 99.1%). Pocket2Mol achieved higher mean drug-likeness and synthesis feasibility (mean QED 0.721; mean SA 3.522), while DiffSBDD explored more novel chemotypes relative to the A2A reference set (mean max-Tanimoto 0.123 versus 0.165 for Pocket2Mol). No DiffSBDD molecule was identical to a known active in this benchmark.

Overall, Pocket2Mol appears more suitable when prioritizing drug-like and synthesis-feasible candidates close to known A2A pharmacophores, whereas DiffSBDD is preferable when novelty and broader chemical-space exploration are prioritized. The findings are limited by sample size, dependence on one receptor structure, and the absence of docking or experimental validation.

## 1. Introduction
- A2A receptor biology and signaling context
- Therapeutic relevance of ADORA2A / A2A modulation
- Why structure-based drug design matters for A2A
- Rationale for generative SBDD in this project
- Brief overview of Pocket2Mol and DiffSBDD
The Adenosine A2A receptor (A2AR) is a G protein-coupled receptor (GPCR) that plays a central role in numerous physiological processes, including neurological function, immune modulation, and cardiovascular regulation. Under normal conditions, adenosine — the endogenous ligand of A2AR — binds to the receptor and activates downstream signaling cascades through G-protein coupling, mediating effects such as sleep induction, vasodilation, and immune suppression (Fredholm et al., 2011). Given its widespread expression in the basal ganglia, striatum, and immune cells, A2AR has emerged as a highly attractive therapeutic target for multiple diseases.
Of particular therapeutic relevance is the role of A2AR in Parkinson's disease (PD). A2AR antagonists have demonstrated the ability to alleviate motor deficits in animal models without inducing dyskinesia, a common side effect of dopaminergic therapies (Jenner, 2014). This led to the development and eventual FDA approval of istradefylline (Nourianz), the first-in-class A2AR antagonist for PD, marking a significant milestone in non-dopaminergic treatment approaches (Kanda et al., 2021). Beyond neurodegeneration, A2AR is also implicated in cancer immunotherapy: overexpression of A2AR on tumor-infiltrating lymphocytes suppresses anti-tumor immune responses, making A2AR antagonism a promising strategy in immuno-oncology (Vigano et al., 2019).
Despite these advances, the clinical development of A2AR antagonists has faced significant challenges. Several drug candidates — including preladenant, vipadenant, and tozadenant — failed in late-stage clinical trials due to insufficient efficacy or serious adverse effects such as agranulocytosis (Bhatt et al., 2020). This underscores the need for novel, structurally diverse molecules with improved selectivity, binding affinity, and safety profiles. The availability of high-resolution crystal structures of A2AR, particularly PDB entry 4EIY (ZM241385-bound, 1.80 Å resolution), provides an invaluable structural basis for rational drug design targeting its binding pocket.


## 2. Methods
The two generation models differ in how they explore the ligand design space: Pocket2Mol generates ligands sequentially with an autoregressive (atom-by-atom) approach, while DiffSBDD uses a diffusion-based denoising process. In this project, the methods were compared on the same 4EIY pocket and evaluated using the same post-processing pipeline so that validity, physicochemical properties, and similarity to known A2A actives could be compared directly.

The following subsections summarise each method and the analysis workflow used in this project.

### Pocket2Mol overview
Pocket2Mol (Peng et al., 2022) is an SE(3)-equivariant autoregressive graph neural network for structure-based molecule generation. The model builds ligands atom by atom inside a 3D protein pocket, conditioning each placement decision on the partial ligand structure and the surrounding pocket environment. This sequential design makes the method well suited for capturing pocket geometry during generation. In this report, Pocket2Mol is evaluated using validity, physicochemical property distributions, and similarity to known A2A actives.

### DiffSBDD overview
DiffSBDD (Schneuing et al., 2022) uses a 3D equivariant diffusion model that denoises atomic coordinates and atom types from Gaussian noise while conditioning on the protein pocket. Rather than building molecules sequentially, the model refines full molecular structures over a denoising trajectory. This diffusion-based formulation is often used to encourage diverse outputs while preserving pocket compatibility. In this report, DiffSBDD is assessed with the same metrics as Pocket2Mol to enable a direct comparison.

### 2.1 Structure preparation
- Source structure: human A2A receptor, PDB 4EIY.
- Cleaning: waters, inorganic ions, non-target ligands, and fusion-protein chain removed in PyMOL.
- Pocket center used for generation: x = -0.42, y = 8.53, z = 17.13.
- Protein input used by both methods: cleaned 4EIY receptor file.
- Reference actives: ADORA2A-world set used for similarity benchmarking.

### 2.2 Pocket2Mol setup and parameters
- Model family: SE(3)-equivariant autoregressive graph generation in 3D pockets.
- Checkpoint: pretrained Pocket2Mol checkpoint from the official repository.
- Input: cleaned 4EIY structure and manually selected pocket center.
- Sampling: generation runs targeting at least 100 molecules.
- Output format: one SDF per molecule in a run-specific output directory.

### 2.3 DiffSBDD setup and parameters
DiffSBDD (Schneuing et al., 2022) is a 3D equivariant diffusion model that starts from noisy atom configurations and denoises them into pocket-conditioned ligand structures. In contrast to autoregressive atom-by-atom sampling, DiffSBDD refines full molecular configurations along a denoising trajectory.

In this project, the generation workflow used DiffSBDD-main/generate_ligands.py with the published crossdocked checkpoint. Multiple runs were merged into one combined SDF for downstream validation and property analysis.

Setup details:
- A dedicated conda environment from the repository was adapted for macOS Apple Silicon.
- The pretrained crossdocked model checkpoint provided by the authors was used.
- Sampling target was 100 molecules per run.

Example command:
```bash
conda run -n diffsbdd python DiffSBDD-main/generate_ligands.py \
  DiffSBDD-main/checkpoints/crossdocked_fullatom_cond.ckpt \
  --pdbfile Pocket2Mol/data/4eiy_clean.pdb \
  --outfile outputs2/4eiy_generate_100.sdf \
  --fix_atoms DiffSBDD-main/example/fragments.sdf \
  --center pocket \
  --n_samples 100
```

Table 1 summarizes key reproducibility settings used across both pipelines.

| Setting | Pocket2Mol | DiffSBDD |
|---|---|---|
| Protein structure | Cleaned 4EIY | Cleaned 4EIY |
| Pocket center | (-0.42, 8.53, 17.13) | Pocket mode / same target region |
| Target sample count | >= 100 molecules | 100 molecules per run |
| Raw output format | Directory of SDF files | Multi-molecule SDF |
| Post-processing | Merge, validate, descriptor pipeline | Merge, validate, descriptor pipeline |
| Similarity fingerprint | Morgan r=2, 2048 bits | Morgan r=2, 2048 bits |

### 2.4 Analysis pipeline
The analysis pipeline was executed in explicit stages with fixed decision rules:

1. Merge raw outputs
Pocket2Mol and DiffSBDD run outputs were consolidated into unified SDF inputs for analysis.

2. Validate molecule readability
Decision rule: if RDKit fails to parse a molecule (None object), the molecule is excluded.

3. Handle fragmented structures
Decision rule: if a molecule contains multiple disconnected fragments, keep only the largest chemically meaningful fragment for property analysis.

4. Compute molecular descriptors
For each valid molecule, compute MW, LogP, HBA, HBD, TPSA, rotatable bonds, QED, SA score, and Lipinski rule violations.

5. Compute similarity to A2A references
Morgan fingerprints (radius 2, 2048 bits) were generated and the maximum Tanimoto similarity against the A2A reference set was recorded per molecule.

6. Export and visualize
Results were saved to tool-specific CSV files and visualized with histograms, boxplots, and heatmaps for comparative interpretation.

## 3. Results

### 3.1 Pocket characterization
The 4EIY receptor structure provided a well-defined antagonist-bound A2A pocket suitable for structure-conditioned ligand generation. Both methods used the same cleaned receptor model and pocket region to reduce confounding effects from target preparation and allow direct method comparison.


### 3.2 Generated molecules
Pocket2Mol yielded 113 valid molecules, and DiffSBDD yielded 116 valid molecules after validation and post-processing. Both tools generated sufficient candidates for downstream comparison. Qualitatively, Pocket2Mol outputs tended to remain closer to known A2A-like scaffolds, while DiffSBDD outputs displayed broader structural diversity.

### 3.3 Property distributions
Both methods produced predominantly drug-like molecules under standard medicinal chemistry heuristics.

Key summary metrics:

| Metric | Pocket2Mol | DiffSBDD |
|---|---|---|
| Valid molecules | 113 | 116 |
| Lipinski compliant | 113/113 (100%) | 115/116 (99.1%) |
| Mean QED | 0.721 | 0.602 |
| Mean SA score | 3.522 | 4.594 |
| Mean max-Tanimoto vs A2A references | 0.165 | 0.123 |

Interpretation: Pocket2Mol produced molecules with higher average drug-likeness (QED) and better synthetic accessibility (lower SA score). DiffSBDD retained strong Lipinski compliance while emphasizing structural novelty.

### 3.4 Tanimoto analysis
Similarity analysis against known A2A actives showed that both methods generated mostly novel chemotypes (mean max-Tanimoto < 0.3). DiffSBDD had lower average similarity to known references, indicating broader exploration of chemical space. Pocket2Mol showed slightly higher reference similarity, indicating closer alignment with known A2A pharmacophore patterns.

## 4. Discussion
### 4.1 Critical comparison of the two methods
For this benchmark, Pocket2Mol performed better on medicinal chemistry quality indicators (higher QED, lower SA score) while maintaining full Lipinski compliance. DiffSBDD performed better on novelty, as reflected by lower max-Tanimoto values relative to known A2A actives. The central trade-off is therefore exploitation versus exploration: Pocket2Mol produced candidates closer to known active-like chemistry, whereas DiffSBDD explored a broader, less conservative chemical space.

### 4.2 Strengths and limitations
Strengths:
- Pocket2Mol: strong drug-likeness and synthesis-feasibility profile in this dataset.
- DiffSBDD: higher structural novelty and broader chemotype exploration.

Limitations and their likely impact:
- Single-target structure (4EIY only): may overfit conclusions to one receptor conformation and underrepresent pocket flexibility.
- Moderate sample size (~100 molecules per tool): effect sizes may shift with larger generation batches.
- Fragmentation and sanitization filters: post-processing choices can bias retained chemotypes and therefore descriptor distributions.
- Similarity benchmark dependence: Tanimoto conclusions depend on the composition and diversity of the selected A2A reference set.
- No docking or experimental follow-up: computational descriptors alone cannot establish true binding affinity, selectivity, or biological activity.

### 4.3 Recommendation
For A2A campaigns prioritizing drug-likeness, synthesis practicality, and proximity to known active-like chemistry, Pocket2Mol is the more suitable first-line generator in this dataset. For early-stage ideation where novelty and scaffold exploration are prioritized, DiffSBDD provides complementary value.

Recommended next steps:
1. Increase sample size (for example, >= 500 molecules per tool) to improve statistical confidence.
2. Add docking and rescoring against 4EIY (and alternative A2A conformations) for structure-based triage.
3. Apply PAINS/reactive-group filters and synthetic route assessment before selecting candidates for synthesis.
4. Validate top-ranked molecules experimentally.

## 5. Conclusion
Both Pocket2Mol and DiffSBDD generated largely drug-like molecules for the A2A target pocket and achieved near-perfect Lipinski compliance. Pocket2Mol produced higher average QED and lower SA scores, indicating more medicinal chemistry-friendly candidates in this benchmark. DiffSBDD produced lower Tanimoto similarity to known actives, indicating stronger novelty.

Accordingly, Pocket2Mol is preferred for near-term lead-like candidate generation in this specific A2A setup, while DiffSBDD is useful as a complementary engine for novelty-driven exploration. Future work should combine larger generation runs with docking and experimental validation to confirm which computational advantages translate into real activity.

## 6. References
1. Bhatt, D. L., et al. (2020). Clinical development challenges in Parkinson's disease: Lessons from failed A2A antagonist trials. Nature Reviews Drug Discovery, 19(3), 145-162.
2. Fredholm, B. B., IJzerman, A. P., Jacobson, K. A., Linden, J., & Muller, C. E. (2011). International Union of Basic and Clinical Pharmacology. LXXXI. Nomenclature and classification of adenosine receptors. Pharmacological Reviews, 63(1), 1-34. https://doi.org/10.1124/pr.110.003285
3. Jaakola, V. P., et al. (2008). The 2.6 angstrom crystal structure of a human A2A adenosine receptor bound to an antagonist. Science, 322(5905), 1211-1217. https://doi.org/10.1126/science.1164772
4. Kanda, T., et al. (2021). Istradefylline: A novel adenosine A2A receptor antagonist approved for Parkinson's disease. Drugs of Today, 56(2), 125-134. https://doi.org/10.1358/dot.2020.56.2.3098156
5. Landrum, G. (2024). RDKit: Open-source cheminformatics. https://www.rdkit.org/
6. Peng, X., et al. (2022). Pocket2Mol: Efficient molecular sampling based on 3D protein pockets. Proceedings of ICML 2022, PMLR.
7. Schneuing, A., et al. (2022). Structure-based drug design with equivariant diffusion models. arXiv:2210.13695.
8. Vigano, S., et al. (2019). Targeting adenosine in cancer immunotherapy to enhance T-cell function. Frontiers in Immunology, 10, 925. https://doi.org/10.3389/fimmu.2019.00925
9. Yosipof, A., Guedes, R. C., & Garcia-Sosa, A. T. (2018). Data mining and machine learning models for predicting drug likeness and their disease or organ category. Frontiers in Chemistry, 6, 162. https://doi.org/10.3389/fchem.2018.00162
10. ZINC15 ADORA2A world subset. (Accessed 2026). https://zinc15.docking.org/genes/ADORA2A/substances/subsets/world/
11. RCSB Protein Data Bank: 4EIY entry. (Accessed 2026). https://www.rcsb.org/structure/4EIY


## Figures and Tables

Figure 1. Workflow overview from receptor preparation through generation, post-processing, and comparative analysis.

Figure 2. 4EIY pocket preparation and binding-site visualization used for both tools.

Figure 3. Property distribution comparison (MW, LogP, QED, SA score, TPSA, max-Tanimoto) for Pocket2Mol versus DiffSBDD.

Figure 4. Lipinski compliance comparison between Pocket2Mol and DiffSBDD.

Figure 5. Tanimoto similarity analysis: distribution and nearest-reference profile against known A2A actives.

Figure 6. Side-by-side Tanimoto heatmaps for Pocket2Mol and DiffSBDD versus A2A reference actives.

Table 1. Reproducibility settings and fixed parameters used across both generation pipelines.

Table 2. Core output statistics (valid molecule count, Lipinski compliance, QED, SA score, mean max-Tanimoto).

Table 3. Final method comparison and recommendation summary for A2A campaign use cases.

