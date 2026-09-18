# Donor allocation and reference-cell budget in MuSiC PBMC simulations

Code, exact cell selections, processed simulation inputs, saved predictions, full-precision summaries and figure source data for a controlled reference-design study. The study separates within-type donor allocation from total reference-cell budget. It includes the original cohort, external three-budget experiment, alternative reference grouping, mechanism diagnostics, arrangement-effect distributions and twelve additional reference-sampling blocks.

Archive DOI: https://doi.org/10.5281/zenodo.22834929

## Download and reproduce

Download `code_results.zip`, `inputs.zip` and `predictions.zip` from the same archive record and extract **all three into one empty directory**. Each ZIP starts at that directory's root; do not add another nested folder. This GitHub repository contains the small source/entry files; the archive supplies the complete code/data snapshot.

Use Python 3.12 with the recorded packages (the original environment was Python 3.12.14). From the extracted directory:

```text
python -m pip install -r reproduction_integrated/environment/python_requirements.txt
python reproduce.py all
```

`python reproduce.py original` recalculates the five original saved-prediction stages and compares their principal MAE/RMSE summaries. `python reproduce.py additional` checks all 645,120 additional prediction rows and 322,560 weighted diagnostics, recalculates every additional-draw summary, runs the separate standard-library summary check and compares all resulting CSVs with the archived values. Fresh result folders preserve the saved outputs.

## Find the results

- Main three-budget results: `results/music_external_review/20260917T133518483857Z/`.
- Original four-allocation gradient: `results/music_gradient_review/20260917T105004747661Z/`.
- Alternative grouping: `results/music_reference_composition_review/20260917T143411436461Z/`.
- Mechanism summaries: `results/music_mechanism_review/20260917T211141303956Z/`.
- Arrangement distributions and Figure S4: `results/music_arrangement_distribution_v1/20260918T153955745453Z/`.
- Additional-draw summaries: `results/music_mc_extension_review/20260918T163322469416Z/`.
- Exact selected cell IDs and complete external inventory: `data/selections/`; exact additional selections and reference/target row maps: `data/processed/music_mc_extension_selection/20260918T153826508686Z/`.
- Figure bundles retain their `*.source_data.csv`, plotting source and scientific labels. Protocols are under `analysis/`; source retrieval and versions are in `reproduction_integrated/environment/`.

`RELEASE_MANIFEST.json` lists every code/result file. `DATA_MANIFEST.json` lists the input and prediction files. [SOURCE_MAP.json](SOURCE_MAP.json) identifies unchanged source hashes and the limited publishing adaptations. [EXPLORATION_HISTORY.md](EXPLORATION_HISTORY.md) preserves earlier unsuccessful and contrary analyses.

## Scope

The tested reproduction route recalculates saved predictions and summaries. It does not refit MuSiC, rebuild the full upstream H5AD, or establish clean-machine/cross-platform fitting equivalence. Original fitting scripts and environment specifications are supplied; their recorded Windows paths require adaptation for a new installation. The additional-draw archive adapter explicitly treats absent installed third-party runtime files as historical provenance, while checking all released scientific inputs, predictions and diagnostics. It preserves the original numerical operations.

All effects describe synthetic mixtures within the fixed donor/target designs. Whole reference-sampling blocks determine conditional Monte Carlo precision; they are not independent biological cohorts. Original three-draw and additional twelve-draw results remain separate, with a labelled combined fifteen-draw summary. A near-zero signed mean does not establish negligible effects in every reference arrangement.

## Attribution and reuse

Author: Yunhao Jiang, School of Biological Sciences, University of Edinburgh. Responsibilities are described in [CONTRIBUTIONS.md](CONTRIBUTIONS.md). See [licenses](LICENSE.md), [source attribution](SOURCES.md), and [CITATION.cff](CITATION.cff). This archive contains the study code and data.
