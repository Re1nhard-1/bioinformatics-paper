source("environment/music/runtime.R")
out <- "results/real_bulk_scale"
stopifnot(!file.exists(file.path(out, "residuals.csv.gz")))
input <- jsonlite::fromJSON("results/music_mc_extension_input/20260918T154429497808Z/input_manifest.json")
cells <- read.delim(file.path(input$data_directory, "cells.tsv"), colClasses = "character")
refs <- read.delim(input$reference_file, colClasses = "character")
ref <- refs[refs$reference_id == "d00_d01_d02_b3_n300_balanced", , drop = FALSE]
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
rownames(cm) <- cm$cell_id
sce <- SingleCellExperiment(assays = list(counts = x), colData = S4Vectors::DataFrame(cm))
bulk <- as.matrix(read.delim(gzfile("data/processed/real_bulk/bulk_counts.tsv.gz"), row.names = 1, check.names = FALSE))
basis <- music_basis(sce, markers = rownames(bulk)[rowMeans(bulk) != 0], clusters = "cell_type", samples = "donor", select.ct = input$cell_types, verbose = FALSE)
saveRDS(basis, file.path(out, "diagnostic_basis.rds"))
common <- intersect(rownames(basis$Disgn.mtx), rownames(bulk))
D <- basis$Disgn.mtx[common, , drop = FALSE]
Y <- sweep(bulk[common, , drop = FALSE], 2, colSums(bulk[common, , drop = FALSE]), "/")
profile <- data.frame(gene_id = common, check.names = FALSE)
for (k in seq_along(input$cell_types)) {
  profile[[paste0("theta_", k - 1L)]] <- basis$M.theta[common, k]
  profile[[paste0("D_", k - 1L)]] <- D[, k]
}
write.csv(profile, file.path(out, "diagnostic_reference_profiles.csv"), row.names = FALSE)
pred <- read.csv("results/real_bulk/predictions_d00_d01_d02.csv", check.names = FALSE)
pred <- pred[pred$reference_id == ref$reference_id, , drop = FALSE]
diag <- read.csv("results/real_bulk/diagnostics_d00_d01_d02.csv", check.names = FALSE)
diag <- diag[diag$reference_id == ref$reference_id, , drop = FALSE]
residuals <- list()
for (i in seq_len(nrow(pred))) {
  row <- pred[i, ]
  qsum <- diag[diag$sample_id == row$sample_id, paste0(if (row$method == "nnls") "nnls" else "weighted", "_coefficient_sum")]
  p <- as.numeric(row[paste0("pred_", 0:5)])
  fitted <- as.vector(D %*% p) * qsum / 100
  y <- Y[, row$sample_id]
  keep <- y > 0
  residuals[[i]] <- data.frame(reference_id = ref$reference_id, sample_id = row$sample_id, method = row$method,
    gene_id = common[keep], observed = y[keep], fitted = fitted[keep], residual = y[keep] - fitted[keep])
}
connection <- gzfile(file.path(out, "residuals.csv.gz"), "wt")
write.csv(do.call(rbind, residuals), connection, row.names = FALSE)
close(connection)
jsonlite::write_json(list(reference_id = ref$reference_id, cell_types = input$cell_types, common_genes = length(common),
  reconstruction = "D * saved proportions * saved coefficient sum / 100; positive bulk genes only; no refitting",
  purpose = "Post-hoc gene-level residual diagnosis after scale adjustments failed"),
  file.path(out, "residual_design.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Saved reference profiles and residuals for", nrow(pred), "existing predictions\n")
