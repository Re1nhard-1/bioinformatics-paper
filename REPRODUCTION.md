# Reproduction routes and script roles

## Supported saved-prediction route

Extract all three ZIP files from this version into one empty directory. Install
`reproduction_integrated/environment/python_requirements.txt`, then run
`python reproduce.py all` from that directory. This recalculates the original
five saved-prediction stages and all additional-draw summary CSVs. It does not
refit MuSiC. Inputs and predictions are unchanged from version 1.0.0.

## Original source scripts

`analysis/prepare_music_gradient.py`, `analysis/launch_music_gradient.py`,
`analysis/review_music_gradient.py` and `analysis/plot_music_gradient.py` are
the original gradient scripts. `analysis/review_music_comparison.py` is the
original endpoint review script. The module `analysis/run_deconvolution_pilot.py`
is included unchanged because these and other scripts import its `sha` and
`save_json` file helpers. Its main function is an earlier normalized-cell pilot;
do not run that pilot to reproduce the paper's raw-count MuSiC results.

For imports/help or plotting/extraction source tools, the larger recorded Python
environment is available via `python -m pip install -r environment/requirements-source-tools.txt`.
This adds matplotlib and h5py to the minimal summary environment. For example,
`python analysis/review_music_comparison.py --help` and
`python analysis/launch_music_gradient.py --help` should display their arguments.

The historical source scripts are not portable archive-level commands. They
expect their original timestamped working layout and source manifests. In
particular, the gradient launcher requires the original endpoint run under
`results/music_comparison/20260917T090859412332Z/` and the recorded R installation;
the preparation script also uses earlier `results/deconvolution_v3_1/` records.
These complete historical working directories and installed R dependencies are
not included. Do not invoke historical fitting scripts directly on this compact
archive and assume that passing a startup check validates a full fit.

The parameterized five-stage runner `reproduction_integrated/code/run.py`
accepts explicit input, Rscript, R-library and official-source paths; inspect
`python reproduction_integrated/code/run.py --help`. Setting up that fitting
environment and repeating full fits on a new machine remain separate steps.
Files below `results/`, including `*.producer.py` and `*.style.py`, are unchanged
execution/provenance snapshots, not standalone entry points; import paths and
root resolution belong to their original execution context. Use scripts in
`analysis/` or the designated reproduction entry instead.

## Version 1.0.1 correction

Version 1.0.0 omitted the gradient launcher and the file-helper module, causing
original source imports to fail. Both are now included byte-for-byte, without
changing scientific calculations, saved predictions, results, figures or
contributions. The repair also specifies the optional h5py dependency and the
distinction between historical sources and the tested saved-prediction route.
`DEPENDENCY_CHECK.json` records source startup and local import checks;
`REPRODUCTION_CHECK.json` records saved-prediction verification. Neither check
claims full upstream data reconstruction or cross-platform MuSiC fitting.
