# Source this file from the project root before running any MuSiC script.
# No user or system configuration is written.
music_project <- normalizePath(getwd(), winslash = "/")
music_library <- file.path(music_project, ".tools/R-library")
if (!dir.exists(music_library)) stop("Run from the bioinformatics-paper project root")
.libPaths(c(music_library, .Library))
stopifnot(length(.libPaths()) == 2L)
suppressPackageStartupMessages(library(MuSiC, lib.loc = music_library))
suppressPackageStartupMessages(library(SingleCellExperiment, lib.loc = music_library))
