# Donor allocation and reference-cell budget in MuSiC

Code and data accompanying **Donor-allocation effects vary with reference-cell budget in MuSiC PBMC simulations** by Yunhao Jiang.

The controlled simulations vary cells contributed by three reference donors while holding their identities and the total cells per type fixed. Larger budgets reduce the average allocation penalty, while effects in individual configurations remain. The measured-bulk assessment did not establish transfer of this trend.

## Contents

- `analysis/`: preparation, fitting, summarization and plotting scripts.
- `results/`: figure source data and selected result tables; `manuscript_tables/` contains the displayed tables, including rounding and explanatory cells.
- `reproduction_integrated/`: parameterized reconstruction and fitting code for the five initial simulation stages.
- `environment/`: recorded Python, R and optional Salmon dependencies.

The complete release is archived at [Zenodo](https://doi.org/10.5281/zenodo.22902701):

- **code_results.zip**: this code, all retained simulation and measured-bulk result tables, design files and figure PDFs.
- **data.zip**: processed count matrices, target truths, reference selections, saved predictions and diagnostic inputs.

Extract both ZIPs into the same empty directory, preserving their relative paths. GitHub provides a smaller browsing copy. Original public datasets and exact versions are listed in [SOURCES.md](SOURCES.md).

## Figures and tables

Install `environment/requirements-deconvolution.lock.txt`. From the extracted directory:

```sh
python analysis/plot_manuscript_figures.py
python analysis/plot_real_bulk.py
python analysis/summarize_salmon_pilot.py --plot-only
```

These commands read saved results for Figures 1–2 and S3–S6; they do not refit MuSiC. Output goes to `figures/` beside the extracted directory. Figure S1/S2 source data and their original plotting scripts are included; the supplied PDFs are the manuscript figures. `results/manuscript_tables/` gives Tables 1 and S1–S18. Full-precision simulation summaries remain in the named `results/music_*_review/` directories. Measured-bulk results are in `results/real_bulk/`, `real_bulk_scale/` and `real_bulk_salmon/`.

## Reproduction

The five initial stages are endpoints, intermediate allocations, thinning, external budgets and alternative grouping. Their saved-prediction summaries can be recalculated without fitting:

```sh
python reproduction_integrated/code/summarize.py --endpoints work/predictions/endpoints --gradient work/predictions/gradient --thinning work/predictions/thinning --external work/predictions/external --alternative work/predictions/alternative --output recomputed
```

For full fitting, `reproduction_integrated/code/run.py --help` accepts explicit input, Rscript, R-library and official MuSiC source paths. Use the pinned MuSiC commit in SOURCES.md and the recorded R dependencies. Additional-draw, diagnostic and measured-bulk source scripts are in `analysis/`; their protocols specify parameters and inputs. Original source scripts retain their study directory conventions. The R scripts expect the project working directory and `.tools/R-library`; the integrated runner accepts an explicit library path. Salmon quantification uses the separate Linux environment and downloads in `environment/salmon/` and `analysis/run_salmon_pilot.sh`.

Saved results were reused for this release. The previously tested five-stage reconstruction and summary route does not establish a new-machine or end-to-end reproduction of every later diagnostic. Raw FASTQs, genomes, installed software, private files, work logs and backups are not included. Earlier normalized-cell exploratory findings are retained in Supplementary Table S3; they are separate from the raw-count MuSiC simulations.

## License and citation

Project code is GPL-3.0-or-later; project results, figures and documentation are CC BY 4.0. Upstream data retain their original terms. See LICENSE.md, SOURCES.md and CITATION.cff.
