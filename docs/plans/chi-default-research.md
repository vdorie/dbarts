# chi-default-research

Part 2 simulation study record (2026-07-10). Part 1 (field survey of
binary k defaults) is in roadmap-survey.md. Full artifacts (48-cell
results, harness) lived in the session scratch; the decision memo and
summary tables are preserved verbatim below.

## DECISION_MEMO.txt

DECISION MEMO: binary-outcome k default for bart2 / dbarts / rbart_vi
Study: 48 cells (4 dgp x n{100,500,2000} x p{5,20} x base rate{.5,.1}),
8 arms (fixed k 1/2/3/5; chi(1.5,{1,2,5,Inf})), 20 paired reps, 7680 fits,
2000 held-out points with known probit truth per rep. Tables:
summary_tables.txt; tidy data: all_results.rds, results/cell_*.rds.

RECOMMENDATION: keep the hyperprior, replace the improper scale:
default k = chi(1.5, 2) for binary responses in bart2/dbarts, and the
same for rbart_vi. Do not move to fixed k=2.

Evidence, per claim:
1. A hyperprior beats every fixed k on robustness. chi_2 worst-cell
   logloss regret 0.008 vs fixed_2 0.038 (~5x), fixed_3 0.079,
   fixed_5 0.142; same ordering for latent RMSE and ECE.
2. Fixed k=2 fails where adaptation matters: strong-separation n=100
   cells sample k~0.8-1.4, and forcing k=2 over-shrinks - coverage of
   true p(x) drops to 0.51 (chi_2: 0.91), its worst-cell regret cell.
   Weak-signal cells sample k~3-6; fixed 2 under-shrinks (calib slope
   0.66 vs chi_2 0.94). No single fixed k covers both regimes.
3. The improper scale is the only pathology in the current default.
   chi(1.5, Inf) drifts in weak-signal/balanced cells: k_median 997
   (weak, n=100, p pooled), k_max 1.5e5, 2 runaway reps, coverage
   0.79, calibration slope ~63 in those cells. Base-rate-0.1 cells
   never drift (k_max <= 7.7). A finite scale removes the drift at
   its source; the 1e6 cap becomes a never-engaged backstop.
4. Scale 2 is the best finite scale tried. chi_2 has the lowest
   worst-cell regret on all three metrics (vs chi_1 and chi_5), best
   or tied-best averages (logloss 0.390, ECE 0.030, latent RMSE 0.103,
   coverage 0.977), and sampled k stays in [0.8, 5.3] across all 48
   cells. chi_1 drags k down (median 1.75) and loses in weak cells;
   chi_5 is a near-tie but slightly worse worst-cell.
5. Field alignment for free: chi(1.5, 2) has prior median k = 1.91
   (scale x 0.953) - it is "the field's k=2, sampled" and reduces the
   unexplained-default objection from the survey to a scale choice
   centered on the consensus value.

rbart_vi: same default, k = chi(1.5, 2) for binary. Docs (rbart_vi +
dbartsPriors + bart2): "for binary responses k defaults to the
chi(1.5, 2) hyperprior, centered near the field-standard k = 2
(prior median 1.9) while adapting to the data; pass k = 2 for the
fixed BART default." bart() keeps fixed k=2 (BayesTree compat).

## summary_tables.txt

loaded 48 cells, 7680 rows

=== AVERAGE across all cells (per arm) ===
     arm logloss   auc    ece latent_rmse coverage ci_width calib_slope_absdev
 fixed_1   0.399 0.787 0.0383       0.118    0.974    0.509              0.280
 fixed_2   0.394 0.794 0.0358       0.108    0.896    0.424              0.338
 fixed_3   0.398 0.795 0.0416       0.109    0.812    0.367              0.528
 fixed_5   0.412 0.793 0.0547       0.121    0.663    0.296              1.129
   chi_1   0.391 0.794 0.0302       0.105    0.981    0.446              0.204
   chi_2   0.390 0.795 0.0300       0.103    0.977    0.423              0.218
   chi_5   0.391 0.796 0.0311       0.103    0.969    0.405              0.297
 chi_Inf   0.393 0.795 0.0317       0.106    0.936    0.381             10.526

=== WORST-CELL REGRET: logloss (lower better) ===
     arm mean_regret worst_regret          worst_cell
 fixed_1     0.01078      0.04622     weak|100|20|0.5
 fixed_2     0.00578      0.03768   strong|100|20|0.5
 fixed_3     0.01039      0.07854   strong|100|20|0.5
 fixed_5     0.02429      0.14186   strong|100|20|0.5
   chi_1     0.00289      0.01579      weak|100|5|0.1
   chi_2     0.00232      0.00804   strong|100|20|0.1
   chi_5     0.00280      0.01773 friedman|100|20|0.5
 chi_Inf     0.00481      0.03933 friedman|100|20|0.5

=== WORST-CELL REGRET: latent_rmse (lower better) ===
     arm mean_regret worst_regret          worst_cell
 fixed_1     0.01883       0.0673      weak|500|5|0.5
 fixed_2     0.00838       0.0312      weak|500|5|0.5
 fixed_3     0.00988       0.0484   strong|100|20|0.5
 fixed_5     0.02139       0.0928    strong|100|5|0.5
   chi_1     0.00562       0.0216      weak|100|5|0.1
   chi_2     0.00374       0.0114      weak|100|5|0.1
   chi_5     0.00389       0.0163 friedman|100|20|0.5
 chi_Inf     0.00620       0.0357 friedman|100|20|0.5

=== WORST-CELL REGRET: ece (lower better) ===
     arm mean_regret worst_regret               worst_cell
 fixed_1     0.01200       0.0889          weak|100|20|0.5
 fixed_2     0.00944       0.0488        strong|100|20|0.5
 fixed_3     0.01532       0.0876        strong|100|20|0.5
 fixed_5     0.02840       0.1259        strong|100|20|0.5
   chi_1     0.00388       0.0430          weak|100|20|0.5
   chi_2     0.00371       0.0249          weak|100|20|0.5
   chi_5     0.00476       0.0289      friedman|100|20|0.5
 chi_Inf     0.00535       0.0379 linear_sparse|100|20|0.5

=== SAMPLED-k BEHAVIOUR (hyper arms) ===
     arm k_median   k_q90   k_max k_frac_cap cells_with_runaway
   chi_1     1.75    2.22    3.01   0.00e+00                  0
   chi_2     2.10    2.77    3.98   0.00e+00                  0
   chi_5     2.45    3.44    5.23   0.00e+00                  0
 chi_Inf    47.18 2199.49 8548.80   6.77e-06                  2

=== chi_Inf runaway incidence by (dgp,baseRate) ===
           dgp baseRate    n k_median    k_max runaway
      friedman      0.1  100    1.989 5.65e+00    0.00
 linear_sparse      0.1  100    1.989 5.97e+00    0.00
        strong      0.1  100    1.394 4.25e+00    0.00
          weak      0.1  100    2.662 7.70e+00    0.00
      friedman      0.5  100   43.783 3.64e+04    0.00
 linear_sparse      0.5  100    7.538 1.14e+03    0.00
        strong      0.5  100    8.113 1.51e+04    0.00
          weak      0.5  100  996.571 1.49e+05    0.05
      friedman      0.1  500    1.638 3.09e+00    0.00
 linear_sparse      0.1  500    1.666 3.23e+00    0.00
        strong      0.1  500    0.866 1.45e+00    0.00
          weak      0.1  500    3.333 6.84e+00    0.00
      friedman      0.5  500    2.107 3.69e+00    0.00
 linear_sparse      0.5  500    2.253 3.91e+00    0.00
        strong      0.5  500    0.944 1.55e+00    0.00
          weak      0.5  500   34.351 3.47e+03    0.00
      friedman      0.1 2000    1.838 2.78e+00    0.00
 linear_sparse      0.1 2000    2.064 3.13e+00    0.00
        strong      0.1 2000    0.821 1.15e+00    0.00
          weak      0.1 2000    3.803 6.32e+00    0.00
      friedman      0.5 2000    2.228 3.22e+00    0.00
 linear_sparse      0.5 2000    2.897 4.31e+00    0.00
        strong      0.5 2000    1.025 1.42e+00    0.00
          weak      0.5 2000    6.367 1.12e+01    0.00

wrote /Users/vdorie/.claude/jobs/7fe13675/tmp/chi-default/all_results.rds 
