# Portable runtime for the unchanged execution-script snapshot.
# Both library and official source are supplied explicitly by code/run.py.
music_library <- Sys.getenv("MUSIC_LIBRARY")
music_source <- Sys.getenv("MUSIC_SOURCE_ROOT")
stopifnot(dir.exists(music_library), dir.exists(music_source))
.libPaths(c(music_library, .Library))
suppressPackageStartupMessages(library(MuSiC, lib.loc = music_library))
suppressPackageStartupMessages(library(SingleCellExperiment, lib.loc = music_library))
stopifnot(as.character(packageVersion("MuSiC")) == "1.0.0")
official <- new.env(parent = asNamespace("MuSiC"))
for (f in c("analysis.R", "construct.R", "utils.R")) {
  sys.source(file.path(music_source, "R", f), envir = official)
}
for (name in c("music_prop", "music_basis", "music.iter", "music.basic")) {
  stopifnot(identical(body(get(name, official)), body(get(name, asNamespace("MuSiC")))))
}
locked <- read.delim("environment/r_packages.tsv", colClasses = "character")
loaded <- installed.packages(lib.loc = .libPaths())
stopifnot(all(locked$Package %in% rownames(loaded)))
stopifnot(identical(unname(loaded[locked$Package, "Version"]), locked$Version))
