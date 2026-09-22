source("environment/music/runtime.R")
prepared <- "data/processed/real_bulk_scale"
source_info <- jsonlite::fromJSON(file.path(prepared, "sources.json"))
input <- jsonlite::fromJSON(source_info$reference_input_manifest)
types <- input$cell_types
cells <- read.delim(file.path(input$data_directory, "cells.tsv"), check.names = FALSE, colClasses = "character")
refs <- read.delim(file.path(prepared, "references.tsv"), check.names = FALSE, colClasses = "character")
genes <- readLines(file.path(input$data_directory, "genes.tsv"))
truth <- read.csv("data/processed/real_bulk/flow_truth.csv", check.names = FALSE)
arms <- c("length_scaled", "author_tpm")
bulks <- lapply(arms, function(arm) as.matrix(read.delim(gzfile(file.path(prepared, paste0(arm, ".tsv.gz"))), row.names = 1, check.names = FALSE)))
names(bulks) <- arms
stopifnot(as.character(packageVersion("MuSiC")) == "1.0.0", all(vapply(bulks, function(x) identical(colnames(x), truth$sample_id), logical(1))))
out <- "results/real_bulk_scale"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
append_csv <- function(x, path) write.table(x, file = path, sep = ",", row.names = FALSE,
  col.names = !file.exists(path), append = file.exists(path), na = "NA")
read_matrix <- function(d) {
  path <- input$donor_matrices$path[input$donor_matrices$donor == d]
  h <- gzfile(path, "rb")
  on.exit(close(h))
  x <- methods::as(Matrix::readMM(h), "CsparseMatrix")
  cm <- cells[cells$donor == d, , drop = FALSE]
  dimnames(x) <- list(genes, cm$cell_id)
  x
}
.music_trace <- new.env(parent = emptyenv())
.music_trace$items <- list()
invisible(trace("music.iter", where = asNamespace("MuSiC"), print = FALSE, exit = quote({
  res <- returnValue()
  st <- get(".music_trace", envir = .GlobalEnv)
  st$items[[length(st$items) + 1L]] <- list(convergence = res$converge, n_features = nrow(D), variance_finite = all(is.finite(res$var.p)))
})))
counts_by_donor <- list()
for (i in seq_len(nrow(refs))) {
  ref <- refs[i, , drop = FALSE]
  prediction_file <- file.path(out, paste0("predictions_", ref$triple_key, ".csv"))
  diagnostic_file <- file.path(out, paste0("diagnostics_", ref$triple_key, ".csv"))
  completed <- character()
  if (file.exists(prediction_file)) {
    previous <- read.csv(prediction_file)
    previous_diag <- read.csv(diagnostic_file)
    stopifnot(all(table(previous$arm) == 24L), all(table(previous_diag$arm) == 12L), setequal(previous$arm, previous_diag$arm))
    completed <- unique(previous$arm)
  }
  if (setequal(completed, arms)) next
  needed <- strsplit(ref$reference_donors, "|", fixed = TRUE)[[1]]
  for (d in needed) if (is.null(counts_by_donor[[d]])) counts_by_donor[[d]] <- read_matrix(d)
  ids <- as.integer(strsplit(ref$reference_columns_R, "|", fixed = TRUE)[[1]])
  cm <- cells[ids, , drop = FALSE]
  x <- do.call(cbind, lapply(needed, function(d) counts_by_donor[[d]][, cm$cell_id[cm$donor == d], drop = FALSE]))
  x <- x[, cm$cell_id, drop = FALSE]
  rownames(cm) <- cm$cell_id
  stopifnot(ncol(x) == 1800L, all(table(cm$donor, cm$cell_type) == 100L), max(abs(Matrix::colSums(x) - as.numeric(cm$total_counts))) == 0)
  sce <- SingleCellExperiment(assays = list(counts = x), colData = S4Vectors::DataFrame(cm))
  for (arm in setdiff(arms, completed)) {
    bulk <- bulks[[arm]]
    .music_trace$items <- list()
    invisible(capture.output(fit <- music_prop(bulk.mtx = bulk, sc.sce = sce, markers = NULL,
      clusters = "cell_type", samples = "donor", select.ct = types, cell_size = NULL,
      ct.cov = FALSE, verbose = FALSE, iter.max = 1000, nu = 0.0001, eps = 0.01,
      centered = FALSE, normalize = FALSE)))
    diag <- do.call(rbind, lapply(.music_trace$items, as.data.frame))
    stopifnot(nrow(diag) == ncol(bulk))
    diag$sample_id <- colnames(bulk)
    diag$reference_id <- ref$reference_id
    diag$arm <- arm
    predictions <- list(weighted = fit$Est.prop.weighted, nnls = fit$Est.prop.allgene)
    rows <- lapply(names(predictions), function(method) {
      p <- predictions[[method]][colnames(bulk), types, drop = FALSE]
      stopifnot(all(is.finite(p)), min(p) >= -1e-10, max(abs(rowSums(p) - 1)) < 1e-10)
      tab <- data.frame(arm = arm, reference_id = ref$reference_id, triple_key = ref$triple_key, block = 3L,
        budget = 300L, level = "balanced", dominant_donor = ref$dominant_donor,
        sample_id = colnames(bulk), method = method, check.names = FALSE)
      for (k in seq_along(types)) tab[[paste0("pred_", k - 1L)]] <- p[, k]
      tab
    })
    append_csv(diag, diagnostic_file)
    append_csv(do.call(rbind, rows), prediction_file)
    cat(i, "/", nrow(refs), ref$reference_id, arm, "complete\n")
    flush.console()
  }
}
