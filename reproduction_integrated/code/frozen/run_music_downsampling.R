# Extension adapter: official package unchanged; paired donor/type thinning, separate protocol.
source("environment/music/runtime.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
input <- jsonlite::fromJSON(args[1])
out <- args[2]
donor <- args[3]
mode <- if (length(args) >= 4L) args[4] else "full"
stopifnot(mode %in% c("full", "replay", "technical"))
technical <- mode == "technical"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
started <- proc.time()["elapsed"]
stopifnot(as.character(packageVersion("MuSiC")) == "1.0.0")
data_dir <- input$data_directory
genes <- readLines(file.path(data_dir, "genes.tsv"))
cells <- read.delim(file.path(data_dir, "cells.tsv"), check.names = FALSE, colClasses = "character")
cases <- read.delim(input$case_file, check.names = FALSE, colClasses = "character")
targets <- read.delim(file.path(data_dir, "targets.tsv"), check.names = FALSE, colClasses = "character")
read_matrix <- function(name) {
  handle <- gzfile(file.path(data_dir, name), "rb")
  on.exit(close(handle))
  methods::as(Matrix::readMM(handle), "CsparseMatrix")
}
counts_all <- read_matrix("counts.mtx.gz")
target_all <- read_matrix("targets.mtx.gz")
stopifnot(nrow(counts_all) == length(genes), ncol(counts_all) == nrow(cells),
          identical(as.integer(cells$matrix_column_R), seq_len(nrow(cells))))
dimnames(counts_all) <- list(genes, cells$cell_id)
dimnames(target_all) <- list(genes, targets$target_name)
cases <- cases[cases$held_out == donor & cases$policy == ifelse(mode == "replay", "replay", "thin5"), , drop = FALSE]
targets <- targets[targets$held_out == donor, , drop = FALSE]
if (technical) {
  cases <- cases[1, , drop = FALSE]
  targets <- targets[1:3, , drop = FALSE]
}
bulk <- as.matrix(target_all[, targets$target_name, drop = FALSE])
types <- input$cell_types
truth <- as.matrix(data.frame(lapply(targets[paste0("true_", 0:5)], as.numeric)))
colnames(truth) <- types
stopifnot(ncol(truth) == 6, max(abs(rowSums(truth) - 1)) < 1e-12)
.music_trace <- new.env(parent = emptyenv())
.music_trace$items <- list()
invisible(trace("music.iter", where = asNamespace("MuSiC"), print = FALSE,
  exit = quote({
    trace_result <- returnValue()
    trace_store <- get(".music_trace", envir = .GlobalEnv)
    trace_store$items[[length(trace_store$items) + 1L]] <- list(
      convergence = trace_result$converge, n_features = nrow(D),
      nnls_coefficient_sum = sum(trace_result$q.nnls),
      weighted_coefficient_sum = sum(trace_result$q.weight),
      variance_finite = all(is.finite(trace_result$var.p)))
  })))
all_checks <- list()
all_diagnostics <- list()
for (j in seq_len(nrow(cases))) {
  cfg <- cases[j, , drop = FALSE]
  ids <- as.integer(strsplit(cfg$reference_columns_R, "|", fixed = TRUE)[[1]])
  ref_meta <- cells[ids, , drop = FALSE]
  expected_cells <- as.integer(cfg$reference_cells)
  expected_per_type <- as.integer(cfg$cells_per_type)
  expected_per_donor_type <- as.integer(cfg$cells_per_donor_type)
  stopifnot(expected_cells %in% c(90L, 360L),
            length(ids) == expected_cells, length(unique(ids)) == expected_cells,
            !(donor %in% ref_meta$donor), length(unique(ref_meta$donor)) == 3L,
            all(table(factor(ref_meta$cell_type, levels = types)) == expected_per_type),
            all(table(ref_meta$donor, ref_meta$cell_type) == expected_per_donor_type))
  ref_counts <- counts_all[, ids, drop = FALSE]
  rownames(ref_meta) <- colnames(ref_counts)
  sce <- SingleCellExperiment(assays = list(counts = ref_counts), colData = S4Vectors::DataFrame(ref_meta))
  .music_trace$items <- list()
  fit <- music_prop(bulk.mtx = bulk, sc.sce = sce, markers = NULL, clusters = "cell_type", samples = "donor",
                    select.ct = types, cell_size = NULL, ct.cov = FALSE, verbose = FALSE,
                    iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE)
  diag_rows <- .music_trace$items
  stopifnot(length(diag_rows) == nrow(targets))
  if (technical) {
    invisible(untrace("music.iter", where = asNamespace("MuSiC")))
    fit_plain <- music_prop(bulk.mtx = bulk, sc.sce = sce, markers = NULL, clusters = "cell_type", samples = "donor",
                       select.ct = types, cell_size = NULL, ct.cov = FALSE, verbose = FALSE,
                       iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE)
    trace_difference <- max(abs(fit$Est.prop.weighted - fit_plain$Est.prop.weighted),
                            abs(fit$Est.prop.allgene - fit_plain$Est.prop.allgene))
    stopifnot(trace_difference == 0)
  }
  methods_to_export <- list(music_weighted = fit$Est.prop.weighted, music_nnls = fit$Est.prop.allgene)
  rows <- list()
  for (method in names(methods_to_export)) {
    pred <- methods_to_export[[method]]
    stopifnot(setequal(colnames(pred), types), setequal(rownames(pred), targets$target_name))
    pred <- pred[targets$target_name, types, drop = FALSE]
    stopifnot(all(is.finite(pred)), min(pred) >= -1e-10, max(abs(rowSums(pred) - 1)) < 1e-10)
    errors <- pred - truth
    frame <- data.frame(case_id = cfg$case_id, held_out = donor, triple_id = as.integer(cfg$triple_id),
      block = as.integer(cfg$block), scenario = cfg$scenario, policy = cfg$policy,
      reference_cells = expected_cells,
      method = method, mixture_id = as.integer(targets$mixture_id), mae = rowMeans(abs(errors)),
      rmse = sqrt(rowMeans(errors^2)), check.names = FALSE)
    for (k in 1:6) {
      frame[[paste0("true_", k - 1L)]] <- truth[, k]
      frame[[paste0("pred_", k - 1L)]] <- pred[, k]
    }
    rows[[method]] <- frame
  }
  write.csv(do.call(rbind, rows), file.path(out, paste0(cfg$case_id, ".csv")), row.names = FALSE)
  diagnostic <- do.call(rbind, lapply(seq_along(diag_rows), function(m) {
    data.frame(case_id = cfg$case_id, mixture_id = as.integer(targets$mixture_id[m]),
               as.data.frame(diag_rows[[m]]), check.names = FALSE)
  }))
  all_diagnostics[[j]] <- diagnostic
  all_checks[[j]] <- list(case_id = cfg$case_id, rows = nrow(targets) * 2L,
      elapsed_seconds = unname(proc.time()["elapsed"] - started), n_common_genes = nrow(fit$Weight.gene),
      maxiter_count = sum(diagnostic$convergence == "Reach Maxiter"))
  if (j %% 10L == 0L || technical) cat(donor, j, "/", nrow(cases), "reference cases; elapsed", round(proc.time()["elapsed"] - started, 1), "s\n")
}
write.csv(do.call(rbind, all_diagnostics), file.path(out, "diagnostics.csv"), row.names = FALSE)
report <- list(status = "completed", donor = donor, technical_only = technical, mode = mode,
   cases = nrow(cases), targets_per_case = nrow(targets), prediction_rows = nrow(cases) * nrow(targets) * 2L,
   all_cases = all_checks, elapsed_seconds = unname(proc.time()["elapsed"] - started))
if (technical) report$traced_vs_plain_maximum_difference <- trace_difference
jsonlite::write_json(report, file.path(out, "fold.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
