source("environment/music/runtime.R")
out <- "results/real_bulk_salmon"
stopifnot(!file.exists(file.path(out, "predictions.csv")))
bulk <- as.matrix(read.delim("data/processed/real_bulk_salmon/bulk_counts.tsv", row.names = 1, check.names = FALSE))
design <- jsonlite::fromJSON("results/real_bulk_scale/residual_design.json")
types <- design$cell_types
basis <- readRDS("results/real_bulk_scale/diagnostic_basis.rds")
common <- intersect(rownames(basis$Disgn.mtx), rownames(bulk))
y <- bulk[common, 1]
y <- y / sum(y)
keep <- y > 0
primary_file <- file.path(out, "saved_basis_fit.rds")
if (file.exists(primary_file)) {
  fit <- readRDS(primary_file)
} else {
fit <- MuSiC::music.iter(y[keep], basis$Disgn.mtx[common[keep], , drop = FALSE],
  colMeans(basis$S, na.rm = TRUE), basis$Sigma[common[keep], , drop = FALSE],
  iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE)
  saveRDS(fit, primary_file)
}
predictions <- list()
collect <- function(arm, weighted, nnls) {
  do.call(rbind, lapply(c("weighted", "nnls"), function(method) {
    p <- if (method == "weighted") as.numeric(weighted) else as.numeric(nnls)
    stopifnot(length(p) == 6L, all(is.finite(p)), min(p) >= -1e-10, abs(sum(p) - 1) < 1e-10)
    row <- data.frame(arm = arm, reference_id = design$reference_id, sample_id = "453W", method = method)
    for (k in 0:5) row[[paste0("pred_", k)]] <- p[k + 1L]
    row
  }))
}
predictions[[1]] <- collect("salmon_saved_basis", fit$p.weight, fit$p.nnls)
diagnostics <- list(data.frame(arm = "salmon_saved_basis", common_genes = length(common),
  positive_genes = sum(keep), convergence = fit$converge, variance_finite = all(is.finite(fit$var.p))))
symbols <- read.delim("data/processed/real_bulk_salmon/tx2gene.tsv")
symbols <- symbols[!duplicated(symbols$gene_id), ]
rpl <- symbols$gene_id[grepl("^RP[LS]", symbols$gene_name)]
quality <- list(saved_basis_shared_genes = length(common), saved_basis_genes = nrow(basis$Disgn.mtx),
  rpl_rps_bulk_fraction_on_primary_genes = sum(y[names(y) %in% rpl]),
  selected_reference = design$reference_id)
input <- jsonlite::fromJSON("results/music_mc_extension_input/20260918T154429497808Z/input_manifest.json")
cells <- read.delim(file.path(input$data_directory, "cells.tsv"), colClasses = "character")
refs <- read.delim(input$reference_file, colClasses = "character")
ref <- refs[refs$reference_id == design$reference_id, , drop = FALSE]
stopifnot(nrow(ref) == 1L)
ids <- as.integer(strsplit(ref$reference_columns_R, "|", fixed = TRUE)[[1]])
cm <- cells[ids, , drop = FALSE]
genes <- readLines(file.path(input$data_directory, "genes.tsv"))
read_part <- function(d) {
  h <- gzfile(input$donor_matrices$path[input$donor_matrices$donor == d], "rb")
  on.exit(close(h))
  x <- methods::as(Matrix::readMM(h), "CsparseMatrix")
  dimnames(x) <- list(genes, cells$cell_id[cells$donor == d])
  x[, cm$cell_id[cm$donor == d], drop = FALSE]
}
x <- do.call(cbind, lapply(unique(cm$donor), read_part))
x <- x[, cm$cell_id, drop = FALSE]
stopifnot(ncol(x) == 1800L, max(abs(Matrix::colSums(x) - as.numeric(cm$total_counts))) == 0)
shared <- intersect(rownames(x), rownames(bulk))
quality$shared_prefilter_genes <- length(shared)
quality$reference_umi_retained_after_prefilter <- sum(x[shared, ]) / sum(x)
rownames(cm) <- cm$cell_id
sce <- SingleCellExperiment(assays = list(counts = x[shared, ]), colData = S4Vectors::DataFrame(cm))
second_basis <- music_basis(sce, non.zero = TRUE, markers = rownames(bulk)[bulk[, 1] > 0],
  clusters = "cell_type", samples = "donor", select.ct = types, cell_size = NULL, ct.cov = FALSE, verbose = FALSE)
second_common <- intersect(rownames(second_basis$Disgn.mtx), rownames(bulk)[bulk[, 1] > 0])
stopifnot(length(second_common) >= 0.2 * min(sum(bulk[, 1] > 0), nrow(sce)))
D2 <- second_basis$Disgn.mtx[second_common, , drop = FALSE]
S2 <- colMeans(second_basis$S, na.rm = TRUE)
Sigma2 <- second_basis$Sigma[second_common, , drop = FALSE]
stopifnot(identical(colnames(D2), types), all(is.finite(D2)), all(is.finite(S2)), all(is.finite(Sigma2)))
y2 <- bulk[second_common, 1]
y2 <- y2 / sum(y2)
second <- MuSiC::music.iter(y2, D2, S2, Sigma2, iter.max = 1000, nu = 0.0001, eps = 0.01,
  centered = FALSE, normalize = FALSE)
saveRDS(second, file.path(out, "shared_prefilter_fit.rds"))
predictions[[2]] <- collect("salmon_shared_prefilter", second$p.weight, second$p.nnls)
diagnostics[[2]] <- data.frame(arm = "salmon_shared_prefilter", common_genes = length(shared),
  positive_genes = length(second_common), convergence = second$converge, variance_finite = all(is.finite(second$var.p)))
write.csv(do.call(rbind, predictions), file.path(out, "predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(out, "diagnostics.csv"), row.names = FALSE)
jsonlite::write_json(quality, file.path(out, "gene_coverage.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
print(do.call(rbind, predictions))
print(do.call(rbind, diagnostics))
