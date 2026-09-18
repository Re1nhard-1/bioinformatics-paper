# MuSiC fixed-budget allocation gradient v1

2026-09-17. Frozen before new intermediate-allocation predictions. Follow-up chosen after observing the original endpoints and thinning experiment; exploratory extension, not preregistration or independent validation. Earlier protocols and results remain immutable.

## Question and scope

With the same three reference donors and 60 cells per type, how do absolute error and the within-method imbalance penalty vary across 20:20:20, 30:15:15, 40:10:10 and 50:5:5? Compare both outputs from unchanged official MuSiC: weighted and internal ordinary NNLS sharing the donor-equal basis. A curve is not assumed monotonic. Do not select favorable donors, targets or methods after fitting.

Reuse all six SLE source donors, six PBMC cell types, all ten reference triples per held-out donor, all three saved reference draw blocks, and all 60 original raw-count synthetic targets per held-out donor. Target compositions are shared across donors, but target cells are donor-specific. Rotate each of the three reference donors as dominant at every unequal allocation. Existing20:20:20/50:5:5 predictions are reused, not refitted except one technical replay. Existing5:5:5 thinning is a separate secondary policy question because it changes total size; exclude it from this fixed-budget curve.

This single extension supersedes the former no-additional-grid investment limit in response to the user's specific request. Do not add other budgets, seeds, algorithms, tissues or unplanned ratios in this run. Earlier literature/source assessments apply; no global novelty claim or inference of clinical sampling thresholds is made.

## Paired construction before fitting

From each original held-out/triple/block group, reconstruct the first50 saved random cells per donor and type using that donor's original dominant50:5:5 case. Check these ordered prefixes against original v3.1 saved cell IDs and against all original shorter subsets. For each dominant donor, take first30/15/15 or40/10/10 cells per type using those same permutations. Each reference contains360 unique cells, 60 per type, with no held-out donor. Donor sets, target bytes, feature input and method settings remain unchanged. A donor-specific prefix is nested as its allotted count changes; the whole reference need not be a subset of the larger-dominance reference because other donor allotments shrink.

There are6 ×10 ×3 ×3 ×2 =1,080 distinct new references, generating129,600 new prediction rows (60 targets ×2 outputs per reference). Reused endpoints contain720 references and86,400 rows. The full four-level analysis has1,800 distinct references and216,000 prediction rows. A balanced reference shared by all three dominant-donor contrasts remains one observation, not three independent replicates.

## Fixed implementation and endpoints

Use official MuSiC1.0.0 commit f21fe67f5670d5e9fca0ad7550abaae3423eb59c, existing locked R environment, flat music_prop, complete original raw-count genes, markers=NULL, iter.max=1000, nu=0.0001, eps=0.01, centered=FALSE, normalize=FALSE, ct.cov=FALSE, cell_size=NULL. Preserve both Est.prop.weighted and Est.prop.allgene. Record effective support, convergence and nonfinite variance; support can change with reference allocation. No outcome-dependent preprocessing or parameter rescue.

Primary descriptive endpoint: MAE in percentage points and MAE(level) minus MAE(balanced), for each method at each of the three unequal levels. Pair at held-out donor/triple/block/target with the same balanced result, average targets, blocks, reference triples and dominant-donor choices equally, then give six held-out donors equal weight. Report the difference between methods' imbalance penalties at every level alongside absolute errors. The original endpoint result remains unchanged.

Report adjacent differences for20→30,30→40,40→50 and whether each is positive, zero or negative, without fitting a chosen curve. Show all six donor profiles. Also report all180 held-out/triple/dominant arrangements averaged over targets and blocks, all18 held-out/block summaries, reversals and ties at each adjacent step, and nondecreasing-profile counts (numerical tie tolerance1e-10 percentage points, not an effect threshold). No selection by sign. Secondary: RMSE, per-type absolute errors, reference-draw SD (ddof=1 over the same three draws, then average). Three draws have limited Monte Carlo precision.

No population P values, paired test, confidence intervals or fitted threshold. Leave-one-donor-out folds share reference donors; 180 configurations and216,000 predictions do not increase biological sample size beyond the six source donors. Descriptive reference-draw SD must not be labelled a population CI. No unplanned effect-size cutoff determines whether to report the results.

## Controls, bounded execution and independent checks

1. Before fitting, verify source/input hashes, all1,080 constructed references, all original prefixes, allocation counts and target identities. Save case selections plus protocol/code hashes; all raw input bytes remain unchanged.
2. Replay the original first balanced reference with all60 targets through the new adapter; both outputs must match original results within1e-10. Run the first30:15:15 case on the first three targets with trace and without trace; predictions must match exactly. These controls cannot be used to change scientific settings.
3. Run the complete grid with at most three R workers and a1,800-second total fitting deadline including controls. Preserve partial/failed outputs. Any package, nonfinite prediction or incomplete-grid failure prevents claiming complete results; do not remove failures or relax tolerances. Iteration-limit events, if any, remain in the report.
4. A separate untraced official caller independently reconstructs all reference IDs and replays one reference at each new level per held-out donor (12 references ×60 targets ×2 outputs =1,440 rows), maximum prediction tolerance1e-10. Independently recalculate error metrics. Root checks all129,600 new saved MAE/RMSE, pairing and aggregate identities. Avoid replaying the full old grid.
5. Produce a figure with all donor profiles, full four-level data, mean curves and paired penalties; line segments connect tested allocations only. Preserve numerical source CSV, PDF/SVG/PNG, caption and code hashes; visually inspect. Update the scientific text whether results are monotonic, flat or reversed. Keep rendered older review snapshots unchanged until later integrated packaging.

Completion: complete verified four-level descriptive comparison and a clear next decision. Then return to independent-data provenance/overlap and budget-availability audit; do not presume more simulated settings resolve external validity or guarantee publication.
