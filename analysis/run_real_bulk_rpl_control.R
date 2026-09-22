source("environment/music/runtime.R")
out <- "results/real_bulk_scale"
stopifnot(!file.exists(file.path(out, "rpl_control_predictions.csv")))
basis <- readRDS(file.path(out, "diagnostic_basis.rds"))
excluded <- readLines(file.path(out, "excluded_rpl_rps_genes.txt"))
design <- jsonlite::fromJSON(file.path(out, "residual_design.json"))
bulk <- as.matrix(read.delim(gzfile("data/processed/real_bulk/bulk_counts.tsv.gz"), row.names = 1, check.names = FALSE))
common <- setdiff(intersect(rownames(basis$Disgn.mtx), rownames(bulk)), excluded)
D <- basis$Disgn.mtx[common, , drop = FALSE]
S <- colMeans(basis$S, na.rm = TRUE)
Sigma <- basis$Sigma[common, , drop = FALSE]
stopifnot(identical(colnames(D), design$cell_types), identical(colnames(D), names(S)))
Y <- sweep(bulk[common, , drop = FALSE], 2, colSums(bulk[common, , drop = FALSE]), "/")
rows <- list()
diagnostics <- list()
for (i in seq_len(ncol(Y))) {
  keep <- Y[, i] > 0
  fit <- MuSiC::music.iter(Y[keep, i], D[keep, , drop = FALSE], S, Sigma[keep, , drop = FALSE],
    iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE)
  outputs <- list(weighted = fit$p.weight, nnls = fit$p.nnls)
  for (method in names(outputs)) {
    p <- outputs[[method]]
    stopifnot(all(is.finite(p)), min(p) >= -1e-10, abs(sum(p) - 1) < 1e-10)
    row <- data.frame(reference_id = design$reference_id, sample_id = colnames(Y)[i], method = method)
    for (k in 0:5) row[[paste0("pred_", k)]] <- p[k + 1]
    rows[[length(rows) + 1L]] <- row
  }
  diagnostics[[i]] <- data.frame(sample_id = colnames(Y)[i], convergence = fit$converge, features = sum(keep), variance_finite = all(is.finite(fit$var.p)))
}
write.csv(do.call(rbind, rows), file.path(out, "rpl_control_predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(out, "rpl_control_diagnostics.csv"), row.names = FALSE)
cat("Completed the fixed RPL/RPS-prefix exclusion control for one reference and 12 participants\n")
