# MuSiC mechanism diagnostics v1

Specified 18 September 2026, after the original/external outcomes and D035 discussion, before any new diagnostic fit. This is a post-hoc descriptive diagnosis, not preregistered confirmation or a change to the frozen scientific endpoint.

## Question and bounded scope

Do larger donor/type quotas show lower observed reference-profile draw dispersion, and how do between-donor Sigma and fitted gene weights vary across the same allocation/budget conditions? These are separate measurements; their association does not establish mediation of MAE.

Use the original external grouping and unchanged official MuSiC commit f21fe67f5670d5e9fca0ad7550abaae3423eb59c, the same 14 donors, six types, three draws, and budgets 60/120/300. Retain only balanced and 10:1:1, including all three dominant-donor choices. There are 504 distinct reference matrices and 2,016 logical held-out configurations. Do not change targets, cell prefixes, genes, method parameters, or donor eligibility.

Select eight existing mixture IDs common to every donor: ID0 (equal composition) plus the seven IDs among1–59 with smallest SHA256 hexadecimal digest of `music-mechanism-v1|mixture|{id}`. Keep all112 selected targets. Selection uses identifiers only, not errors or diagnostic outcomes. Each reference predicts32 targets from its four eligible held-out donors:16,128 weighted fits and the corresponding16,128 internal ordinary outputs. All original840-target results remain authoritative for the main endpoint.

## Correct baseline and metrics

The exported93,892-cell input is a union of reference prefixes and target cells, not all114,524 eligible donor/type cells. No exported-union profile is treated as full-inventory truth. No new source download is required. Source inventory N comes from the frozen cohort_inventory.csv, not the exported counts.

1. For each donor/type at n=5,10,20,25,40,50,100,250, reconstruct each of the three exact nested-prefix profiles theta=sum(counts by gene)/sum(all counts), on all30,172 genes. Record the trace of sample covariance across draws, sum_g var_b(theta_gb), with denominator2; mean pairwise total-variation distance over the three pairs; and CV of the three mean library sizes. Record source N and n/N. Profiles reused across reference triples are deduplicated only after checking equality.
2. For each reference/type record the sum over genes of official between-donor Sigma, mean library size, and squared norm of the official donor-equal mean profile. Reference-all-zero genes contribute zero to theta/Sigma summaries. Independently reconstruct theta, S, M.theta, Sigma and design matrix from raw selected counts and compare with music_basis at tolerance1e-10.
3. For each fitted target/reference record support size, normalized-weight effective fraction `(1/sum(p^2))/G`, and top-ceiling(1% of G) normalized-weight share, where p is weights divided by their sum on that fit's support.
4. For each unequal reference paired to balanced at the same triple/draw/budget/target, record support intersection/union and Jaccard; Spearman correlation of positive finite weights on identical shared gene IDs; total-variation distance between separately normalized weights on that intersection; and the fraction of each fit's total weight retained on the intersection. Missing fitted weights are not zero-filled. We do not reconstruct iteration denominator components from final residuals.

The major sampling-precision summaries are the paired quota trajectories5→10→25 (minor donors),20→40→100 (balanced),50→100→250 (dominant). Report all84 donor/type endpoint ratios and direction counts, plus donor means and type summaries. Three draws give imprecise conditional dispersion estimates; no population CI or P value is attached. Finite-population sampling fractions can approach1 and must be reported. No25-cell threshold is tested.

Sigma allocation contrasts are unequal-minus-balanced within reference/type, with equal type/draw/triple/dominant-choice weighting and then equal donors. Weight-comparison summaries similarly give equal targets, draws, triples and dominant choices within donor, then equal donors. Retain individual records and contrary directions. Do not subtract profile dispersion from Sigma as a purported estimate of biological variance. Do not correlate thousands of overlapping rows and report them as independent observations.

## Technical gate before full diagnostic dispatch

The first lexicographically ordered triple, block0, all three budgets, balanced plus ratio10 with its first lexicographic donor, supplies six reference controls. Compare its32-target call to four separate8-target calls, both output methods and actual gene weights on matched gene IDs/NA masks. Maximum finite difference must be≤1e-10, support masks identical. Compare every saved control prediction with its frozen original full-grid prediction, max≤1e-10. Check basis reconstruction as above. Stop on failure; preserve output and source snapshots; do not weaken criteria after seeing results. A corrected implementation uses a new versioned filename.

After controls pass, run all14 triples with at most3 workers. Shared wall-clock limit30minutes for fitting, output budget1GiB; stop on errors/nonfinite diagnostics/iteration-limit fits or control failure. Preserve incomplete outputs and report them, rather than silently skipping cases. Compare all32,256 new prediction rows to their frozen original counterparts and recompute MAE/RMSE from six proportions. Duplicate profile diagnostics and structural row counts must agree. Source hashes and script snapshots accompany each invocation.

## Evidence and stopping rules

This milestone ends with complete tables, one accurately labelled diagnostic figure, and an evidence review. A reduction in observed draw dispersion supports only that conditional observation. Sigma/weight changes alone do not establish that sampling precision caused the MAE budget contrast. No new real-bulk experiment, alternative salt search, marker selection, or component ablation is included. Component swapping would need a separately justified protocol; no automatic extension to obtain a preferred story.

New diagnostic files do not overwrite old scientific results, protocols, snapshots or figures. Main-paper integration follows review of the complete diagnostic outcomes, including results that do not support the candidate explanation.
