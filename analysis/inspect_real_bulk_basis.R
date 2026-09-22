source("environment/music/runtime.R")
input <- jsonlite::fromJSON("results/music_mc_extension_input/20260918T154429497808Z/input_manifest.json")
cells <- read.delim(file.path(input$data_directory, "cells.tsv"), colClasses = "character")
refs <- read.delim(input$reference_file, colClasses = "character")
ref <- refs[refs$reference_id == "d00_d01_d02_b3_n60_balanced", , drop = FALSE]
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
common <- intersect(rownames(basis$Disgn.mtx), rownames(bulk))
D <- basis$Disgn.mtx[common, , drop = FALSE]
y <- bulk[common, "453W"]
y <- y / sum(y)
keep <- y > 0
q <- stats::coef(nnls::nnls(D[keep, , drop = FALSE], y[keep]))
p <- q / sum(q)
saved <- read.csv("results/real_bulk/predictions_d00_d01_d02.csv", check.names = FALSE)
saved <- saved[saved$reference_id == ref$reference_id & saved$sample_id == "453W" & saved$method == "nnls", ]
difference <- max(abs(p - as.numeric(saved[paste0("pred_", 0:5)])))
markers <- c(CD3D = "ENSG00000167286", CD3E = "ENSG00000198851", MS4A1 = "ENSG00000156738", CD79A = "ENSG00000105369", LYZ = "ENSG00000090382", NKG7 = "ENSG00000105374", ACTB = "ENSG00000075624", B2M = "ENSG00000166710")
tab <- data.frame(gene = names(markers), ensembl = markers, bulk_relative = y[markers], check.names = FALSE)
for (ct in input$cell_types) tab[[ct]] <- basis$M.theta[markers, ct]
write.csv(tab, "results/real_bulk/basis_marker_diagnostic.csv", row.names = FALSE)
all_predictions <- do.call(rbind, lapply(list.files("results/real_bulk", pattern = "^predictions_.*[.]csv$", full.names = TRUE), read.csv))
near_zero <- lapply(split(all_predictions, all_predictions$method), function(g)
  list(B_zero = sum(g$pred_0 <= 1e-10), T_zero = sum(g$pred_2 + g$pred_3 <= 1e-10), fit_count = nrow(g)))
jsonlite::write_json(list(posthoc_reason = "B and T estimates were near zero in the initial real-bulk assessment",
  reference_id = ref$reference_id, sample_id = "453W", gene_alignment_by_exact_identifier = TRUE,
  independent_nnls_max_difference = difference, nnls_cell_types = input$cell_types,
  independent_nnls_proportions = p, reference_cell_size = basis$M.S,
  near_zero_tolerance = 1e-10, near_zero_estimates = near_zero),
  "results/real_bulk/input_diagnostic.json", auto_unbox = TRUE, pretty = TRUE, digits = 16)
stopifnot(difference < 1e-10)
print(tab)
print(p)
