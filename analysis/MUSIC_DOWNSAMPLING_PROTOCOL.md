# MuSiC donor-balanced thinning v1

2026-09-17. Frozen before any new thinning predictions. This is an exploratory follow-up chosen after the original fixed-budget results, not public preregistration or an independent replication. It does not change `MUSIC_COMPARISON_PROTOCOL.md` or replace any earlier run.

## Reason and literature gate

The user asked whether the contribution could be strengthened and then requested the next step after the proposal to assess keeping all available imbalanced cells versus donor-balanced thinning. Submission preparation is paused for this single bounded comparison. The earlier prohibition on arbitrary grid expansion remains; the new policy question is separately versioned and authorized by this continuation.

Huuki-Myers et al. (2025), DOI 10.1186/s13059-025-03552-3, already studied equalizing cell-type proportions by downsampling to the rarest type (Bisque/hspe, real bulk/RNAScope evaluation; Supplementary Figs S29–S30). Its donor subsampling varies donor number. Neither finding makes generic downsampling harm novel. The proposed distinction is thinning donor-by-type groups within a fixed three-donor MuSiC reference. MuSiC already gives donor-level signatures equal weight; thinning is not a correction for unequal donor weights. A literature audit must permit this narrow contrast before any fitting; its full evidence is in `research/DOWNSAMPLING_LITERATURE_AUDIT.md`. No absent exact parameter match will be called proof of novelty.

## Available-data decision

For each original 50:5:5 donor allocation within each of six types, compare:

- **Keep all:** existing official MuSiC outputs using 60 cells/type, 360 reference cells total.
- **Thin to five per donor/type:** retain the first five cells in each saved donor/type random permutation, giving 5:5:5, 15 cells/type and 90 reference cells total. Every cell is an actual subset of the available 50:5:5 reference; no new cells, labels or draws are acquired. This discards 75% of reference cells.

The 50:5:5 inventory is simulated from a larger deposited dataset; it is not an observed clinical resource constraint. The total sample size deliberately differs: that loss is part of the available-data policy, not a confound to conceal. Do not interpret its effect as a pure causal effect of balance. The original 20:20:20 case is an **expanded-inventory comparator only**: it cannot be constructed from the two donor/type groups with five cells in the 50:5:5 inventory. Show it only with that explicit label; do not recommend unavailable cells as preprocessing.

## Fixed data and scope

Reuse exact original counts and metadata from `results/music_input/20260917T090438982999Z/input_manifest.json`, all six held-out donors, all ten reference triples, all three existing reference blocks and all 60 saved raw-count targets per held-out donor. Target compositions are shared across donors; target cells are donor-specific. No new expression processing, feature panel, method, seed, tissue or allocation grid is added. The same exploratory SLE cohort is not external validation.

Metadata-only feasibility `results/downsampling_feasibility/20260917T100531718232Z/audit.json` established that all three dominant-donor choices within each donor/triple/block share the identical first-five subset. Therefore 540 imbalanced references need only **180 distinct thinning fits**: 6 × 10 × 3. With 60 targets and two official outputs there are **21,600 new prediction rows**. Pair each subset with all three parent allocations in summaries; repeated subset values do not create independent evidence.

Use unmodified official MuSiC 1.0.0 commit `f21fe67f5670d5e9fca0ad7550abaae3423eb59c`, the existing locked R environment, flat `music_prop`, full original-gene inputs and `markers=NULL`. Preserve all other defaults from the original protocol (iter.max=1000, nu=0.0001, eps=0.01, centered=FALSE, normalize=FALSE, cell_size=NULL, ct.cov=FALSE). Export both `Est.prop.weighted` (primary) and `Est.prop.allgene` (secondary matched method output). Gene support, library-size and variance estimates may change with thinning; this is a total workflow effect.

## Endpoints and pairing

For each prediction compute six-type MAE and RMSE exactly as in the original study; multiply proportion-scale errors by 100 when reporting percentage points. Primary estimand is **weighted MAE(thin) − weighted MAE(keep all)**. Positive means thinning increases error. Within each held-out donor, average equally over 60 targets, three existing blocks, all ten reference triples and three dominant donors. Give the six held-out donors equal weight in the final mean.

Report both policies' absolute MAE, each donor's paired contrast, all 180 donor/triple/dominant-donor contrasts averaged over blocks and targets, and all 18 donor/block summaries. Report any reversals, ties or failures. Secondary outcomes: analogous ordinary-NNLS error and RMSE; per-type absolute errors; conditional reference-draw SD using ddof=1 across the same three blocks for fixed donor/triple/policy/target/type, then averaging. Three draws have low Monte Carlo precision. No population P values, confidence intervals, effect-selected subgroup or outcome-dependent draw exclusion. No choice of favored policy from synthetic truth will be used to construct the input.

## Controls, independent checks and stopping

1. Verify original input/output/source hashes used by this extension and lock new case selections and script/protocol hashes before launch. Check all first-five cells against the saved original v3.1 donor/type lists, subset membership in each parent, identical reference donor set, held-out exclusion, six types and exact group sizes. Confirm no matrix or target bytes change.
2. First run a technical control: replay the first original balanced reference with its complete 60-target batch, comparing both outputs to saved original results with maximum absolute tolerance 1e-10. This checks the new orchestration independently of new outcomes. Also run the first thinning case on three targets, compare traced and untraced output exactly and check diagnostics. These are controls, not a subset used to choose analysis settings.
3. Run all 180 new references only after controls pass. A single new R script changes input/reference-size assertions and output labels; official MuSiC code and parameters remain unchanged. Use at most three concurrent R processes and a 1,800-second **total** fitting deadline. Preserve partial outputs on failure; never claim full completion if any case fails. Do not relax numerical checks to rescue a result.
4. Independently check selected inputs/targets and at least one thinning case per held-out donor with all 60 targets through a separately written untraced R caller, tolerance 1e-10. Recalculate all saved MAE/RMSE and summary contrasts. Preserve all convergence/support diagnostics, including iteration-limit or nonfinite cases.

Completion is a complete descriptive policy comparison and an explicit contribution decision, regardless of sign. No publication success is defined by a favorable result, significance, or an arbitrary one-percentage-point cutoff. Assess whether the observed magnitudes, heterogeneity and prior-work overlap justify a separately specified independent-data check. Do not automatically add a grid or rewrite the old manuscript as a stronger validated claim. The current Word/PDF review snapshot remains unchanged until the new evidence is evaluated.
