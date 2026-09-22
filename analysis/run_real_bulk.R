source("environment/music/runtime.R")
args <- commandArgs(trailingOnly = TRUE)
worker <- if (length(args)) as.integer(args[1]) else 0L
workers <- if (length(args) > 1) as.integer(args[2]) else 1L
stopifnot(worker >= 0L, worker < workers, as.character(packageVersion("MuSiC")) == "1.0.0")
prepared <- "data/processed/real_bulk"
sources <- jsonlite::fromJSON(file.path(prepared, "sources.json"))
input <- jsonlite::fromJSON(sources$reference_input_manifest)
types <- input$cell_types
cells <- read.delim(file.path(input$data_directory, "cells.tsv"), check.names = FALSE, colClasses = "character")
refs <- read.delim(input$reference_file, check.names = FALSE, colClasses = "character")
genes <- readLines(file.path(input$data_directory, "genes.tsv"))
bulk <- as.matrix(read.delim(gzfile(file.path(prepared, "bulk_counts.tsv.gz")), row.names = 1, check.names = FALSE))
truth <- read.csv(file.path(prepared, "flow_truth.csv"), check.names = FALSE)
stopifnot(identical(colnames(bulk), truth$sample_id), ncol(bulk) == 12L)
out <- "results/real_bulk"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_matrix <- function(path) {
  h <- gzfile(path, "rb")
  on.exit(close(h))
  methods::as(Matrix::readMM(h), "CsparseMatrix")
}
append_csv <- function(x, path) write.table(x, file = path, sep = ",", row.names = FALSE,
  col.names = !file.exists(path), append = file.exists(path), na = "NA")
.music_trace <- new.env(parent = emptyenv())
.music_trace$items <- list()
invisible(trace("music.iter", where = asNamespace("MuSiC"), print = FALSE, exit = quote({
  res <- returnValue()
  st <- get(".music_trace", envir = .GlobalEnv)
  st$items[[length(st$items) + 1L]] <- list(convergence = res$converge, n_features = nrow(D),
    variance_finite = all(is.finite(res$var.p)), nnls_coefficient_sum = sum(res$q.nnls),
    weighted_coefficient_sum = sum(res$q.weight))
})))
triples <- sort(unique(refs$triple_key))
triples <- triples[(seq_along(triples) - 1L) %% workers == worker]
started <- proc.time()["elapsed"]
for (triple in triples) {
  selected <- refs[refs$triple_key == triple, , drop = FALSE]
  stopifnot(nrow(selected) == 96L)
  prediction_file <- file.path(out, paste0("predictions_", triple, ".csv"))
  diagnostic_file <- file.path(out, paste0("diagnostics_", triple, ".csv"))
  completed <- character()
  if (file.exists(prediction_file)) {
    previous <- read.csv(prediction_file, check.names = FALSE)
    previous_diag <- read.csv(diagnostic_file, check.names = FALSE)
    stopifnot(all(table(previous$reference_id) == 24L), all(table(previous_diag$reference_id) == 12L),
      setequal(previous$reference_id, previous_diag$reference_id))
    completed <- unique(previous$reference_id)
  }
  if (length(completed) == nrow(selected)) next
  needed <- strsplit(selected$reference_donors[1], "|", fixed = TRUE)[[1]]
  counts_by_donor <- lapply(needed, function(d) {
    path <- input$donor_matrices$path[input$donor_matrices$donor == d]
    stopifnot(length(path) == 1L)
    x <- read_matrix(path)
    cm <- cells[cells$donor == d, , drop = FALSE]
    stopifnot(nrow(x) == length(genes), ncol(x) == nrow(cm))
    dimnames(x) <- list(genes, cm$cell_id)
    x
  })
  names(counts_by_donor) <- needed
  for (i in seq_len(nrow(selected))) {
    ref <- selected[i, , drop = FALSE]
    if (ref$reference_id %in% completed) next
    ids <- as.integer(strsplit(ref$reference_columns_R, "|", fixed = TRUE)[[1]])
    cm <- cells[ids, , drop = FALSE]
    quota <- table(cm$donor, cm$cell_type)
    expected <- if (ref$level == "balanced") rep(as.integer(ref$minor_count), nrow(quota)) else
      ifelse(rownames(quota) == ref$dominant_donor, as.integer(ref$major_count), as.integer(ref$minor_count))
    stopifnot(length(ids) == as.integer(ref$budget) * 6L, !anyDuplicated(ids),
      setequal(unique(cm$donor), needed), all(quota == expected), all(colSums(quota) == as.integer(ref$budget)))
    x <- do.call(cbind, lapply(needed, function(d) counts_by_donor[[d]][, cm$cell_id[cm$donor == d], drop = FALSE]))
    x <- x[, cm$cell_id, drop = FALSE]
    rownames(cm) <- cm$cell_id
    stopifnot(max(abs(Matrix::colSums(x) - as.numeric(cm$total_counts))) == 0)
    sce <- SingleCellExperiment(assays = list(counts = x), colData = S4Vectors::DataFrame(cm))
    .music_trace$items <- list()
    invisible(capture.output(fit <- music_prop(bulk.mtx = bulk, sc.sce = sce, markers = NULL,
      clusters = "cell_type", samples = "donor", select.ct = types, cell_size = NULL,
      ct.cov = FALSE, verbose = FALSE, iter.max = 1000, nu = 0.0001, eps = 0.01,
      centered = FALSE, normalize = FALSE)))
    diag <- do.call(rbind, lapply(.music_trace$items, as.data.frame))
    stopifnot(nrow(diag) == ncol(bulk))
    diag$sample_id <- colnames(bulk)
    diag$reference_id <- ref$reference_id
    predictions <- list(weighted = fit$Est.prop.weighted, nnls = fit$Est.prop.allgene)
    rows <- lapply(names(predictions), function(method) {
      p <- predictions[[method]][colnames(bulk), types, drop = FALSE]
      stopifnot(all(is.finite(p)), min(p) >= -1e-10, max(abs(rowSums(p) - 1)) < 1e-10)
      tab <- data.frame(reference_id = ref$reference_id, triple_key = triple, block = as.integer(ref$block),
        budget = as.integer(ref$budget), level = ref$level, dominant_donor = ref$dominant_donor,
        sample_id = colnames(bulk), method = method, check.names = FALSE)
      for (k in seq_along(types)) tab[[paste0("pred_", k - 1L)]] <- p[, k]
      tab
    })
    append_csv(diag, diagnostic_file)
    append_csv(do.call(rbind, rows), prediction_file)
    if (i %% 8L == 0L) {
      cat(triple, i, "/", nrow(selected), "references; worker", worker, "; elapsed", round(proc.time()["elapsed"] - started), "s\n")
      flush.console()
    }
  }
  rm(counts_by_donor, x, sce)
  invisible(gc())
}
cat("Worker", worker, "complete\n")
