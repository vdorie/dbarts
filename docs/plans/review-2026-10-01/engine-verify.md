# Review 3 - engine lens - independent verification

Tree 01dee4b4 (pinned worktree, read only). Library r3-lib for the shipped behaviour; an instrumented
copy of the same tree (git archive, scratch only) in r3-verify-engine-lib with a cone census and a
runtime switch R3_FIX that replaces the cone score with a log-space reference. All probes are mine,
written from scratch; scratch dirs r3-verify-engine-{cone,onecut,census,tn,src}.

## engine-01 - CONFIRMED, widened - MAJOR

(a) Score error. Verbatim copy of monotoneAdaptiveSimpson, monotoneIntegrate and coneProbability
compiled against R's Rf_pnorm5 (r3-verify-engine-cone/cone.cpp), unbounded pair with unit joint sd,
log(engine) - pnorm(-gap, log.p = TRUE):

    sR/sL:    0.05   0.2    1     1.2    5     20     50
    gap 6     0.000  0.000  0.000 0.000  0.000  0.000  0.000
    gap 7     0.001  0.019 -0.077 0.077 -0.331  0.021 -0.287
    gap 7.5   0.003  0.012  0.034 0.181 -1.042 -3.868 -3.868
    gap 10    0.009 -0.032 -0.189 -0.005 -1.388 -4.275 -4.275
    gap 30    0.010  0.030 -0.180 0.055  0.504 -14.13 -23.69
    gap 37   -0.019 -0.047 -0.190 -0.174 -1.994 -2.372 -34.88
    gap 40    -Inf everywhere (exact about -804)

First gap with |error| > 1e-6 is 5.0-5.8 for every ratio, > 1e-3 at 5.6-6.6, > 0.1 nats at 6.5-7.1 once
sR/sL >= 1. Cause as claimed: the absolute tolerance (1e-12/16 per panel) is met without refinement
once the whole integral is below ~1e-12, so the 16-panel estimate stands; its error scales with how
narrow the peak is (width sL / sqrt(sL^2 + sR^2) in the u units against panels 4.75 wide), so
sR/sL < 1 stays within ~0.05 nats until the underflow. The ~38 sd underflow is exp underflow
(peak height e^(-gap^2 / 2)), not only the +-38 clamp.

NEW second mechanism (not in engine.md): when the lower leaf has a finite lower bound aL from a frozen
neighbour, coneProbability forms the inner mass as gaussianCdf((min(bL, x) - mL) / sL) - lowerL in
lower-tail arithmetic. Once (aL - mL) / sL exceeds ~8.3 both CDFs round to 1 and the integrand is 0
everywhere, so the integral is exactly 0 and twoLeafCoupledLogMarginal returns the -HUGE_VAL
"infeasible" sentinel for a feasible move with a log mass near -42, not -745. Captured from a real
fit ("dip" below): lowR = aL = -0.054755, bR = bL = Inf, mL = -0.129222, sL = 0.008561,
mR = -0.059441, sR = 0.019915: engine 0; log-space reference -41.875; R integrate in upper-tail
arithmetic -41.876.

(b) Sampler. My own one-cut probe (r3-verify-engine-onecut/probe.R): one binary constrained
predictor, one tree, fixed sigma, deterministic data, n0 rows at x = 0 and n1 at x = 1 sitting delta
below. The tree prior ratio was not assumed: the unconstrained twin and a mild-gap constrained
control pick it (CGM with both children owing 1 - 0.95/4; the twin 0.9153 vs 0.9115, the controls
leaf 0.7714 vs 0.7706, joint 0.6312 vs 0.6268, gap 2.4 leaf 0.5528 vs 0.5508). 20000 draws:

    n0 n1 delta sigma  gap  sU/sL | leaf: exact / engine-law / sampler | joint: exact / engine-law / sampler
    200 1  4.2  0.5   8.29  13.5 | 0.167 / 0.542 / 0.545             | 0.091 / 0.372 / 0.365
    200 1  3.8  0.5   7.48  13.4 | 0.197 / 0.004 / 0.004             | 0.109 / 0.002 / 0.002
    200 1  3.1  0.2  15.40  14.0 | 0.059 / 0.306 / 0.302             | 0.030 / 0.181 / 0.182
     50 1  1.7  0.2   8.27   6.9 | 0.167 / 0.430 / 0.434             | 0.091 / 0.274 / 0.278
    200 5  1.9  0.5   8.39   6.2 | 0.121 / 0.322 / 0.322             | 0.065 / 0.192 / 0.193

The sampler follows the quadrature's law, not the exact one, and the error goes either way. With
R3_FIX the same sampler gives 0.169 / 0.096 (row 1) and 0.202 / 0.110 (row 2): exact.

Why moderate split probabilities survive at gaps of 8-15 sd: the free split's likelihood gain is
about gap^2 / 2 and the cone mass costs about the same, so the move is decided by the O(1) remainder,
exactly where the 0.1-35 nat errors land. Small cone masses are not "decided anyway".

(c) Both priors are affected identically (same score; "leaf" only adds the order-count ratio).

(d) Real fits, census over every cone evaluation (n = 500, x1 constrained increasing, 3 predictors,
noise sd 0.5, 500 + 500 sweeps; "joint" 75 trees / "leaf" 20 trees):

    scenario                      cone calls  mass<1e-12   |err|>0.01  >1 nat  false zero  mean|err|
    increasing truth                  19443        0             0        0        0       0
    flat in x1                        15528        0             0        0        0       0
    5 outliers of 8 sd                18359        0             0        0        0       0
    steep (10 x1)                     24236        4             4        1        0       0.0001
    decreasing truth (joint)          14178     6648          6534       10        0       0.061
    decreasing truth (leaf)            3765     1974          1917        2        0       0.069
    dip (local decrease, joint)       17963     1505          1469        5        2       0.117
    decreasing truth, n = 2000        13938    11334         10754      129        3       0.344

The unbounded case (closed form) is 45-99% of cone calls (81-99% at 75 trees). The one-leaf normalMass path never erred
(0 of ~200k). Posterior consequence, 6 seeds x {shipped, R3_FIX}, 2000 kept draws, n = 500: no
detectable change in x1 split counts (z 0.8, -0.5, -0.3) or sigma; fitted values move 0.03-0.06
posterior sd on average, 0.10-0.25 at most, within what 6-seed t-statistics allow.
n = 2000 (where 81% of cone calls are below 1e-12): again no detectable change (x1 splits z -1.0 and
-0.5, sigma z 0.2 and 0.7; fitted values 0.06 / 0.13 posterior sd on average, |t| > 3 on 33 and 30 of
2000 rows, about the ~27 that 6-seed t-statistics give under no change).

So: well-specified, flat and outlier fits do not reach the regime; fits whose truth runs against the
constraint do, on half to four fifths of their cone evaluations, with errors that are mostly small in
those designs but large (and sign-random) whenever the touched children differ much in size.

(e) Fix direction: correct, and the log-space part is necessary, a relative tolerance alone is not
sufficient (it repairs neither the lower-tail cancellation nor the exp underflow). Validated: the
R3_FIX log-space reference reproduces the exact one-cut law above. Sketch for an implementer:

1. Replace coneProbability with logConeProbability returning the log mass; twoLeafCoupledLogMarginal
   adds it and returns -HUGE_VAL only on an empty domain (lowR > bR), never on a computed zero.
2. All four bounds infinite (45-99% of calls): Rf_pnorm5((mR - mL) / sqrt(sL^2 + sR^2), 0, 1, 1, 1).
   This also takes most of the cost in TODO monotone-leaf-quadrature off the profile.
3. Otherwise integrate in u (upper leaf standardized) logf(u) = -u^2/2 - log sqrt(2 pi) +
   logStandardNormalMass(zA, (min(bL, mR + sR u) - mL) / sL), zA = (aL - mL) / sL or -Inf, over
   [max((lowR - mR) / sR, (aL - mR) / sR), (hiR - mR) / sR], no +-38 clamp. logf is log-concave
   (a normal log density plus a log normal mass of a concave nondecreasing argument), so factor the
   peak-relative machinery out of monotoneInvertLogConcave (golden-section mode, cut at 50 nats,
   16 panels a side, scale min(1, sL / sR)) into a shared monotoneLogIntegrateLogConcave that returns
   peak + log(integral of exp(logf - peak)), tolerance relative to that integral (O(peak width)).
4. Precondition: logStandardNormalMass returns NaN when hi is one ulp above lo, because Rf_pnorm5 is
   not ulp-monotone: d = farTail - nearTail comes out slightly positive and log(-expm1(d)) is NaN
   (175 of ~1e6 near-equal pairs in [-8, 8]; it hit my reference on a real fit's integration
   endpoint). Clamp d = std::min(d, 0.0). Today it is reachable only through
   drawPairUpperByInversion's density at u within an ulp above aL (rare; MINOR on its own), but
   step 3 evaluates it at exactly that endpoint.
5. oneLeafLogMarginal: base + logStandardNormalMass((a - m) / s, (b - m) / s) (robustness only).

Tests that would have caught it:
- tests/cpp test_monotone.cpp: logConeProbability against the closed form over gap in {0, 2, 5,
  6.5, 7, 7.5, 8, 10, 15, 30, 37, 40, 60, 200} x sR/sL in {0.05, 0.2, 1, 5, 20, 50}, |error| <=
  1e-9 max(1, |exact|); bounded cases with (aL - mL) / sL in {6, 8.7, 20} and finite bL, lowR, hiR
  against a brute-force log-space Simpson (1e5 nodes) in the test; driven through
  logLikelihoodForBranchWithParams on a hand-built two-leaf tree so the seam, not only the helper,
  is covered. Plus logStandardNormalMass(lo, nextafter(lo, Inf)) never NaN for 1e5 lo in [-8, 8].
- Exact gate (monotone-exact-enumeration.R) extension: contrary-data one-cut designs, e.g. 200 rows
  at x = 0, 1 row at x = 1 sitting 4.2 below, sigma 0.5 fixed (gap 8.3, sU/sL 13.5; exact P(split)
  leaf 0.167, joint 0.091; shipped engine 0.545 / 0.365, ~4 s per prior at 20000 draws, fits quick
  mode); one at sU/sL near 1 and gap 12 (the r < 1 control); and a two-cut design where a frozen
  neighbour bounds the lower child ~9 sd above its posterior (the false-zero path). Its reference
  orderProbability must move to log space (closed form for one cut, log-space propagation
  otherwise) or it returns -Inf in the same regime and cannot see the defect.
- SBC cannot see it (data drawn from the monotone prior); a misspecified-truth SBC-style arm is
  not a calibration check, so the exact gate is the right home.

Severity: MAJOR, not BLOCKER. The documented claim (each prior targeted exactly) is false in a
reachable regime and a wrong-sign structure posterior is easy to build with a lone contrary point
next to a large cell; but in the realistic misspecified fits measured here the fitted function did
not move detectably.

## KNOWN truncated-normal-upper-tail - CONFIRMED (new evidence holds) - MINOR

R emulation of ext_rng_simulateTruncatedNormalScale1's bulk branch, mean 0, 2e5 draws, bias in
units of the truncated sd: (7.5, 8.0] 282 values, -0.009; (8.0, 9.0] 7 values, +0.33, KS D 0.14;
(8.2, 9.0] 2 values (8.21 and the clamp at 9), +3.3, KS D 0.42; from a = 8.3 the gap is 0 and Robert
rejection takes over (exact). The defect window is lower bounds about 7.7-8.3 sd above the mean on an
interior category (OrdinalResponse::drawLatents, and the bridge's AL::ordinal at
R_interface_bartcore.cpp ~7055); the top category's one-sided draw is unaffected. So the record's
"precision" is a discretization bias, but only in that window. Fix as TODO says: reflect when
lower > mean (draw -X on (-upper, -lower] about -mean), or build the bulk branch from upper-tail
pnorm/qnorm when a > 0; the existing equivalence re-record applies. Test: KS of 1e5 draws on
(8.0, 9.0] and (8.2, 9.0] against the exact truncated law.
