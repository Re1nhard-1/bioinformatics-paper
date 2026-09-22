# External PBMC budget/allocation study v1

2026-09-17. Frozen before external MuSiC predictions. This is an outcome-unseen external follow-up designed after the original SLE results and external metadata audit, not a preregistered study or proof of person-level independence. Preserve all earlier protocols/results.

## Question, cohort and implementation

Does reference-donor allocation sensitivity persist, diminish or reverse as the total reference budget increases on a fixed new cohort? Use the Perez/CELLxGENE public processed release c55dc602-d168-4d15-acc1-5de4f2f5d551, with the audited raw.X and raw.var, not normalized X. The source also links healthy accession GSE137029; do not say the H5AD was downloaded as a GEO GSE174188 count attachment.

Include all 14 normal/na donor labels meeting at least 250 cells in each of B, cM, T4, T8, ncM, NK after the eight old-source labels were screened out. Inventory: results/cellxgene_inventory/20260917T113340800025Z. All have one sample and one suspension UUID; some span processing cohorts and every library pools donors. Retain those covariates. Do not split libraries into biological subjects, pool different disease states, relabel cells, or rank donors by expression or prediction. The high-inventory inclusion condition and the normal-versus-original-SLE setting limit generalization; matching-label exclusions cannot rule out aliases.

Map the six author labels to the original six analysis types in the original ordering: B cells, CD14+ Monocytes, CD4 T cells, CD8 T cells, FCGR3A+ Monocytes, NK cells. The labels represent corresponding broad populations; do not assert identical annotation truth across studies. Supply all 30,172 unique Ensembl IDs, preserve raw counts and full library sizes, and do not replace IDs with duplicated symbols or intersect the genes with the old study.

Use unchanged official MuSiC 1.0.0 commit f21fe67f5670d5e9fca0ad7550abaae3423eb59c and the existing R lock: flat music_prop, markers=NULL, select.ct=the six types, cell_size=NULL, ct.cov=FALSE, iter.max=1000, nu=0.0001, eps=0.01, centered=FALSE, normalize=FALSE. Save both Est.prop.weighted and Est.prop.allgene (ordinary NNLS sharing the MuSiC basis), effective support and convergence. No new algorithm, marker panel, preprocessing rescue or outcome-based exclusion.

## Bounded reference design

Order the 14 donor labels by ascending SHA256 of UTF-8 `music-external-v1|ring|{donor}`. For the donor at ring position i as held-out target, list the other 13 donors in forward cyclic order i+1,...,i+13. Omit the final donor (i+13) from reference construction for this held-out fold, and partition the first 12 into four consecutive triples. This gives four triples per held-out donor; each available reference donor is used once in that fold except the predeclared omitted donor. Across all 14 folds every donor is omitted once, contributes to references in 12 other folds, and contributes to targets once. Save the full ring and incidence tables; no search for favorable triples. This is a restricted reference-composition design, not all possible triples.

Use three reference-cell draws (blocks 0,1,2). For each donor, type and block, independently permute its sorted source-row indices and keep the first 250. The same prefix is reused wherever that donor/type/block appears in any reference, and across all budgets/allocations. No per-reference cell replacement. This intentionally shares draws across folds/triples and does not create independent replicates.

Budgets B=60,120,300 mean total cells per type across the three reference donors; whole references contain 6B=360,720,1800 cells. At each budget use:

| Allocation level | B=60 | B=120 | B=300 |
| --- | --- | --- | --- |
| balanced | 20:20:20 | 40:40:40 | 100:100:100 |
| 2:1:1 | 30:15:15 | 60:30:30 | 150:75:75 |
| 4:1:1 | 40:10:10 | 80:20:20 | 200:50:50 |
| 10:1:1 | 50:5:5 | 100:10:10 | 250:25:25 |

Rotate all three donors as dominant at every unequal level. A balanced reference is constructed once for each triple/block/budget and shared by its three contrasts. Use each donor/type prefix to its required length. Budget increases at a given allocation create nested cell sets; changes in allocation shrink some donor subsets and grow others, so complete references need not be nested across imbalance levels.

There are 14 × 4 × 3 × 3 × (1+3×3) = 5,040 logical held-out reference cases, each evaluated against 60 targets and two outputs, for 604,800 planned scientific prediction rows. The cyclic construction reuses 14 distinct donor triples; the same 1,260 triple/block/budget/allocation reference matrices each apply to four held-out donors. A later launcher may fit those 240-target batches to avoid rebuilding identical bases, provided every target belongs to an eligible held-out donor and results are exported to the same 5,040 logical cases. The scientific comparisons and counts cannot change through batching.

## Synthetic targets and deterministic sampling

Reuse only the 60 saved 300-cell composition count vectors from results/music_input/20260917T090438982999Z/input_manifest.json and its targets.tsv (all six original donors have the same truth). Do not reuse original cell identities or expression. There is one equal-composition target and 59 original random composition vectors. Validate exact integer counts summing to 300 before constructing external targets. Observed component maxima are 179,159,201,185,201,200 in the six-type order, all below the candidate inventory minimum 250.

For each external donor and target 0,...,59, sample without replacement within each type to its prescribed count, using new donor-specific draws. Cells can recur between different synthetic targets. Sum their full raw count vectors into a gene-by-target count matrix; truth is the six integer counts/300, not RNA fraction. All budgets, triples, dominant choices and reference blocks use exactly the same saved 60 targets for a held-out donor. Donors serving as targets never appear in the same fit's reference; cells may have target and reference roles in different leave-one-donor-out fits. Record this dependence rather than claiming independent folds.

Sampling algorithm: NumPy Generator(PCG64(seed)), where seed is the unsigned big-endian integer of the first 16 bytes of SHA256 of UTF-8 key. Reference key: `music-external-v1|reference|{donor}|{author_type}|{block}`; target key: `music-external-v1|target|{donor}|{author_type}|{target_index}`. Starting pools are source H5AD row indices sorted ascending. Reference uses permutation(pool)[:250]; targets use permutation(pool)[:required_count]. Save row selections and source cell IDs. Do not search or retry seeds based on results.

## Endpoints and interpretation

All endpoints are descriptive. For method m and budget B, the strongest imbalance penalty I_m(B) is MAE(10:1:1) minus MAE(balanced), in percentage points, averaged equally over targets, blocks, four triples and their three dominant donors, then over 14 held-out donors. Primary budget contrast: I_weighted(300) minus I_weighted(60); negative values mean a smaller penalty at larger budget in this design. Report the full 14 paired donor differences and signs regardless of direction. No minimum effect cutoff controls reporting.

Required accompanying results: absolute MAE for both methods at all three budgets/four levels; each level's imbalance penalty and weighted-minus-NNLS penalty difference; corresponding 120-versus-60 and 300-versus-120 contrasts; all donor profiles, each dominant-donor contrast and allocation reversals. A sensitivity penalty must never be reported without absolute error, since better baseline accuracy can coexist with a larger penalty. Secondary: RMSE, type-specific absolute error, reference-draw SD across the same three blocks (ddof=1), clearly conditional with only three draws. No fitted universal budget threshold, population confidence interval, ordinary paired test or P value based on shared folds/targets/reference cells. Computation rows do not increase the biological sample size.

Do not append a new thinning arm or simultaneously expand disease groups. External targets remain synthetic single-cell sums, not measured bulk samples with experimental truth. Cross-study changes in disease/processing/annotation/gene universe are not controlled causal contrasts. Integrate null/reversed outcomes and revisited limitations into the eventual manuscript.

## Input verification, technical gate and execution bounds

1. Before any external predictions, hash this protocol, code and original metadata/selection inputs; save all 5,040 cases, ring/incidence, donor covariates, reference prefixes and 840 target cell lists. Independently verify allocation quotas, nested prefixes, target counts, complete pairing and donor exclusion. No matrix-performance filter.
2. Retrieve the union of selected cells from fixed-version raw.X using identity-checked ranges and existing cache; coalescing contiguous byte ranges is allowed. Validate every selected raw value as finite, nonnegative and integer, every index in [0,30172), CSR row/column sizes, no repeated gene index within a cell, positive library counts, original row-to-cell/metadata mapping and unique IDs. Sort sparse indices if needed without changing counts; never silently merge duplicated gene coordinates. Cross-check a fixed subset of rows directly against HDF5 and verify pseudobulk reconstruction and exact truth counts. Keep Ensembl IDs and raw library totals.
3. Technical gate only in the input-preparation milestone: replay the old first balanced reference with its first three saved targets through unchanged official arguments, tolerance1e-10. On external first ring donor/first triple/block0, fit first three targets for B60 balanced, B120 4:1:1 with first reference donor dominant, and B300 10:1:1 with the same dominant donor. These are 18 unique external prediction rows, recorded as technical controls rather than a completed scientific grid. Trace/plain official calls must agree exactly. Report failures and diagnostics; do not use control error direction to alter the protocol.
4. Input preparation: disk check, less than 8GiB additional network data, at most three concurrent fetches, 45-minute acquisition deadline; failures retain cache and require a documented operational continuation. Do not fabricate a whole-H5AD hash from cached ranges. Do not load normalized X.
5. Scientific full run is a later milestone after the gate: maximum three R workers and 10,800 seconds fitting deadline. Preserve failed/partial outputs and iteration-limit events. No parameter relaxation or omission of failed cases. Any operational revision gets a separate record; any scientific design change requires a new protocol version before the affected predictions.

Completion of this milestone is frozen design, verified complete selected counts/targets and passed technical controls, not external biological validation or a manuscript acceptance claim.
