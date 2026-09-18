# Exploratory history and contrary evidence


The MuSiC comparison followed earlier observations and is an exploratory follow-up. Its protocol was frozen before empirical MuSiC predictions, but the study was not publicly preregistered. Table S3 records why earlier observations cannot be silently repackaged as confirmatory evidence for the current claim. Pancreas data came from four GSE84133 human donors [S4]; the old analyses used five eligible pancreas types, normalized-cell targets and selected marker panels. Earlier marker ranking used additional cells from the training donors beyond the 60-cell signature budget; that earlier budget was not a complete acquisition budget. Their estimates are not pooled with the six-type PBMC raw-count comparison.

**Table S3. Completed exploration and consequences for interpretation.** “Equal” and “pooled” in the first four rows refer to earlier reference averaging rules, not the two official MuSiC outputs. E1–E7 identify exact evidence entries below.

| Stage | Executed scope | Result retained | Consequence and evidence |
| --- | --- | --- | --- |
| Pancreas pilot | Four donors; 30,240 predictions; normalized NNLS | Strong-allocation MAE was 5.786 pooled versus 5.438 equal; 4/12 dominant-donor arrangements worsened. SD was a post hoc diagnostic. | No universal accuracy benefit; E1. |
| Pancreas v2 | Four panels, two solvers, 24 draws; 645,120 predictions | Top100 NNLS equal-minus-pooled MAE was −0.285 points, but removing five abundant genes reversed it to +0.182. The reversal occurred in all three seed blocks. Across eight settings, conditional SD differences were +0.993 to +1.907. | Preserve feature sensitivity and reversal; E2. |
| External PBMC v3 failure and v3.1 | Six donors; two panels/two solvers; full v3.1 grid of 2,419,200 predictions | Initial SLSQP oracle discrepancy 6.722 × 10⁻⁵ exceeded 5 × 10⁻⁵. Uniformly tightening `ftol` from 10⁻¹² to 10⁻¹⁴ yielded a complete passing run; scientific design and acceptance limits stayed fixed. Top100 NNLS MAE improved 0.433 points while conditional SD increased 1.694. | Failed partial outputs retained; the pancreas deletion reversal was not reproduced across the four PBMC settings; E3–E4. |
| Finite-population calibration | Reconstructed saved v2 reference signatures; no new proportion predictions | Strong-allocation equal/pooled signature-variance trace ratio was 3.466 theoretically and 3.411 empirically. | Consistent with known sampling arithmetic, not a new discovery or an explanation of a percentage of proportion SD; E5. |
| Official MuSiC follow-up | Same saved PBMC cell IDs; raw counts; 720 references; 86,400 predictions | Weighted MAE 5.634→6.361; ordinary MAE 9.110→9.227; difference in differences +0.610 points. Weighted absolute error and conditional SD remained lower. | Current narrow additional control; no new algorithm, mechanism, universal threshold or method ranking; E6–E7. |

In the calibration, independent sampling within each donor–type inventory gives Var(mean) = (1/n − 1/N)S². For reference weights w(d), the signature-coordinate variance is Σ_d w(d)²(1/n(d) − 1/N(d))S(d)². Summing coordinates gives a covariance trace without requiring independence between genes. This concerns normalized reference signatures, not the variance of nonlinear estimated proportions. The initial “equal-donor instability is a new finding” framing was abandoned; no negative result or failed run was removed to support the present narrower report.

Historical record identifiers (these identify the original working records; this compact release does not include all earlier exploratory outputs or reports):

- **E1:** `results/deconvolution_pilot/20260916T211211342333Z/run.json`; `research/DECONVOLUTION_PILOT_REPORT.md`.
- **E2:** `results/deconvolution_v2/20260916T213110522685Z/run.json`; `results/deconvolution_v2_review/20260916T213405002442Z/review.json`; `research/DECONVOLUTION_V2_REPORT.md`.
- **E3:** `results/deconvolution_v3/20260916T220825964629Z/failed.json`; `analysis/DECONVOLUTION_V3_1_PROTOCOL.md`; `results/v3_solver_diagnostic/20260916T221515861510Z/run.json`.
- **E4:** `results/deconvolution_v3_1/20260916T221536546175Z/run.json`; `results/deconvolution_v3_review/20260916T222031951122Z/review.json`; `results/v3_independent_review/20260916T221714283707Z/verification.json`.
- **E5:** `results/reference_variance_calibration/20260916T220104708874Z/run.json`; `research/REFERENCE_VARIANCE_CALIBRATION.md`.
- **E6:** `analysis/MUSIC_COMPARISON_PROTOCOL.md`; `results/music_input/20260917T090438982999Z/input_manifest.json`; `results/music_comparison/20260917T090859412332Z/run.json`.
- **E7:** `results/music_review/20260917T091456882105Z/review.json`; `results/music_independent_review/20260917T091408058654Z/final_verification.json`.
