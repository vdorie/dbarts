# Ensemble-scale SBC for every shipped family (2026-08-24, tip b102e17c)

Method. `git archive HEAD` staged to a scratch tree, `R CMD INSTALL --preclean`
into a private lib, every R call prefixed `R_LIBS=<lib>`. The machinery is
benchmarks/R/sbc.R UNCHANGED (source()d, so its CLI block never fires); the only
re-parameterisation is `nTrees`, through a 90-line driver (sbc-logs/run-arm.R).
Scale: n.trees = 200, `bart()`'s shipped default (`dbartsControl` defaults to
75, so 200 brackets both), against the 50 every recorded SBC ran at; n = 150
(160 grouped, 200 bcf), p = 3. Convention from sbc-calibration.md and
sbc-family-tiers.md, unchanged: R = 200 replications, L = 200 retained draws
(150 for the four tier arms), thin = 30, each family's own measured burn in
absolute sweeps (`sbcBurnSweeps`; 72000 for the runSbc arms). Admission, per
sbc-family-tiers.md step 4: the ecdf-diff band Bonferroni'd over the matrix's
functional count - 0.05/83 = 6.02e-4, band ~0.137. PASS/NOTE/FLAG = inside 5%
band / over 5% but inside Bonferroni / over it. All harness self-checks passed
at 200 trees, and `sbc.R discrete-selfcheck` replays BYTE-IDENTICAL.

Wall clock 3 h 43 m (16:11-19:54), 9.45 h of CPU at 3 concurrent; the driver
and band recomputation are in ./sbc-logs/ (raw logs and rank matrices
recording not retained).

## Family table (verdicts at the 5% band; 3 of 83 FLAG at the Bonferroni band)

    family      covered  scale         worst functional  ecdf/band05  verdict  s     log
    gaussian    yes      m=200 n=150   f.star5 0.1028    1.12         6/7 P    3035  gaussian.log
    student(t)  yes      m=200 n=150   avg.f 0.0725      0.78         4/4 P     569  t.log
    + weights   yes      m=200 n=150   avg.f 0.0657      0.72         7/7 P    2500  weighted.log
    probit      yes      m=200 n=150   avg.f 0.0885      0.97         6/6 P    3190  probit.log
    logistic    yes      m=200 n=150   f.star2 0.0717    0.78         6/6 P    3657  logistic.log
    ordinal     yes      m=200 n=150   eta.star3 0.1065  1.15         9/10 P   1445  ordinal.log
    nbinom      yes      m=200 n=150   r 0.2000          2.16         1/3 P    3860  nbinom.log
    multinom    yes      m=200x3 n=150 f.1.1 0.0779      0.84         6/6 P    1262  multinom.log
    bcf         yes      m=200x2 n=200 prog3 0.1261      1.38         12/15 P  3434  bcf.log
    grouped     yes      m=200 n=160   sigma 0.2916      3.18         9/10 P   2579  grouped-gaussian.log
    grouped/pro yes      m=200 n=160   f.star5 0.0822    0.90         9/9 P    2686  grouped-probit.log
    aft         NO       -             -                 -            -        -     see below
    hazard      NO*      -             inherited         -            -        -     see below
    hurdle      NO*      -             inherited         -            -        -     see below
    hetero      NO       -             -                 -            -        -     see below

    ladder points (A4e protocol, thin 90 = 3x spacing, R = 100)
    gaussian    f.star5   1.12 -> 0.47   SHRINKS   1524  gaussian-thin90.log
    ordinal     eta.star3 1.15 -> 0.67   SHRINKS   2155  ordinal-thin90.log
    bcf         prog3     1.38 -> 0.51   SHRINKS   2115  bcf-thin90.log
    bcf         prog1     1.11 -> 1.04   FLAT             (same run)

## The flags and their adjudication

FLAG 1 - nbinom `r` and `agg.psi`. ecdf-diff 0.2000 and 0.1785 against a 0.137
Bonferroni band, chisq p 0.000 both, `avg.mu` (the identified mean) PASSES at
0.0884. `r`'s histogram is a U with a 41/10 spike in the TOP rank bin (mean rank
88.3 of a target 75), `agg.psi` mirrors it with a 44/10 spike in the BOTTOM bin
(mean 64.9): theta0's r sits above nearly all posterior draws and its psi below,
i.e. posterior r too SMALL and psi too LARGE, with mu = r exp(psi) exact.
ADJUDICATION: already discharged, not re-laddered. This reproduces the recorded
flag at four times the forest size. The 2026-08-18 SBC-deepening pass took it to
a full-R thin=300 point (spacing-invariant) and then to an independent
derivation that found the r conditional correct term by term and DEMONSTRATED
the mechanism: at n = 150 the r-vs-psi ridge is effectively non-ergodic
(integrated autocorrelation 2199-6359 sweeps; two chains held disjoint r basins
for 100000 sweeps with identical avg.mu), the cold start at the grid median
making the error one-sided exactly as seen here; the n = 20 control passes
everything. A KNOWN MIXING LIMITATION at ensemble n, not a draw-law defect; new
here is that it is forest-size-invariant, as a ridge the trees never touch is.

FLAG 2 - grouped-gaussian `sigma`. ecdf-diff 0.2916, chisq p 0.000, mean rank
135.2 of 100, histogram monotone increasing (8 -> 41): theta0's sigma above the
posterior, i.e. posterior sigma too small. ADJUDICATION: HARNESS ARTIFACT,
confirmed here rather than assumed. A grouped sampler refuses setResponse, so
the fit is REBUILT per replication and the gaussian response-adaptive scaling
re-anchors to range(y0) while theta0's f0 came from the fixed-scale generator;
the width mismatch lands entirely in the residual variance. Three confirmations
in this run: (a) the same arm's tau (0.0873), b1 (0.0362), b2 (0.0703) and all
five f* PASS - the grouped-tau target itself calibrates; (b) grouped-PROBIT,
which has no response scaling, is 9/9 PASS with tau at 0.0542; (c) plain
gaussian, which reuses one sampler at a fixed scale, has sigma at 0.0513.
Recorded at m = 50 as the same artifact (0.213). Not a sampler finding; it is
an argument for the scale-pin below.

NOTE 3 - bcf `prog1` (0.1014) and `prog3` (0.1261), `b1.minus.b0` (0.1261). All
three are inside the Bonferroni band; all five `eff` cells, `abs.a`, `abs.diff`
and `sigma` PASS. `b1.minus.b0`'s histogram is a U (29 bottom, 21 top): the
documented sign-symmetry - a*mu and b_z*tau are invariant under a joint sign
flip under a symmetric Cauchy a prior, so the RAW glue posteriors are bimodal
and their SBC is ILL-POSED BY CONSTRUCTION, which is why the harness also
reports the identified |a| and |b1-b0| (both PASS). The prog cells carry a mild
high-rank tilt (mean 110.1, 111.7 of 100). Laddered at thin 90: prog3 shrinks
1.38 -> 0.51 and b1.minus.b0 1.38 -> 0.85, prog1 stays flat at 1.11 -> 1.04 with
a healthy chisq p 0.371 there. Reading: the recorded (a, mu-amplitude) ridge
(A4e's H-MIX, bcf-ridge-interweaving) still grazes the prognostic function at
one of five evaluation points. NEW AND GOOD:
the BCF `sigma` channel - the one that read as a calibration DEFECT at m = 50
before A4e (ecdf 0.1344 at thin 120) - is CLEAN at m = 200 at both spacings
(0.0609, 0.0925). The ridge's grip is weaker at the shipped default forest size.

NOTE 4 - gaussian `f.star5` (0.1028, chisq p 0.204) and ordinal `eta.star3`
(0.1065). Both inside the Bonferroni band, both SHRINK under 3x spacing (to 0.47
and 0.67 of band), and the identity of the marginal cell does not survive a
re-run: ordinal's recorded m = 50 flag was gamma3, 0.0739 PASS here, and at the
ladder point the marginal ordinal cell moves again to p4 (1.06) while eta.star3
falls. Five NOTEs over 83 functionals is the 4.15 expected at 5%. Noise.

NO PERSISTENT FLAG THAT IS A DRAW-LAW DEFECT CANDIDATE. Both hard flags are
pre-adjudicated (one a demonstrated mixing limitation with the derivation
checked, one a harness scaling artifact reproduced three ways), and every NOTE
either shrinks under spacing or is ill-posed by construction. At four times the
forest size, gaussian, student-t, weighted gaussian, probit, logistic, ordinal,
multinomial (INCLUDING the raw f_ik cells), grouped random intercepts under both
bases, and BCF's identified effect and prognostic functions all calibrate.

One recorded OPEN item changes status. multinomial's three raw per-forest
`f_ik` cells - the finding sbc-family-tiers.md could not discharge (persistent
U, chisq p 0.000, 0.86-0.91 of band at m = 50, NOT shrinking at 3x spacing,
pre-registered suspect `MultinomialForestCombiner::afterCombine`'s
level-centering precision) - are CLEAN at m = 200: 0.0779 / 0.0701 / 0.0439
against a 0.0924 5% band, chisq p 0.326 / 0.395 / 0.958, no U. The only change
is the forest size, and the per-tree leaf scale falls as 1/sqrt(m), so an
approximated centering precision would show exactly this m-dependence. It does
not clear the approximation - the exact-conditional derivation that plan asked
for remains the right instrument - but it bounds its reach: at the SHIPPED
default the raw levels calibrate.

## Not covered, with the cost to cover

aft - NOT COVERED, and the sharpest gap: AFTResponse (model.hpp:3794) owns real
sampling code (a truncated-normal latent redraw per censored row per sweep), so
it reduces to nothing. The block: censoring STATUS is fixed at creation
(setResponse refreshes bounds but cannot change which rows are censored,
model.hpp:3868) while the generative model makes status_i = 1{logT*_i <= c_i} a
function of theta0. Two enablers: (A) a status setter - engine + bridge + R5
method + refusal parity + tests, ~100-150 lines over three layers, and NEW
PUBLIC SURFACE, so it lands before 1.0-0 or not at all; (B) a creation-time
leaf-scale pin plus rebuild-per-replication - the FIX-A shape
runsbcbcf-repair.md priced at ~30-60 engine lines and did not take (FIX-B
landed). B is cheaper, adds no R surface, and would ALSO retire the grouped
sigma artifact, which is the same re-anchoring. Harness either way ~80 lines. A
zero-censoring AFT arm is free and worthless: with no censored rows the fit is
bit-identical to a gaussian fit on log T, so the gaussian arm IS it.

hazard - NOT COVERED DIRECTLY, CALIBRATION-INHERITED. The recorded reason stands
(the person-period design, hence the tree prior's cut grid, depends on y0). But
hazard adds no engine code - R/dbarts.R:487-536 expands and remaps to the binary
token before the sampler exists - and hazard-reduction.R gates DRAW-FOR-DRAW
bitwise equality with a hand-expanded bart2(probit)/bart2(logistic) fit at 40
trees. The probit and logistic arms above plus that reduction ARE ensemble-scale
evidence for hazard. A direct arm needs the at-risk set frozen (no censoring,
all K periods observed), ~60 harness lines; low value given the reduction.

hurdle - NOT COVERED DIRECTLY, CALIBRATION-INHERITED for the two components.
bart2Hurdle (R/bart.R:2206) literally calls bart2(probit) and bart2(gaussian),
and hurdle-reduction.R gates them draw-for-draw at 30 trees, so the probit and
gaussian arms cover both components at ensemble scale. NOT covered: the
combine/retransform step, whose analytic oracle is still owed (hurdle.md step
2). A direct arm hits the same y0-dependent-design block (the y > 0 subset); the
owed oracle, not SBC, is the right instrument.

hetero - NOT COVERED, the other gap with real sampling code (the multiplicative
variance roll, chain.hpp:724) and no reduction argument. Prior draws never reach
varianceForest_: there is no sampleTreesFromPrior for the variance forest, so
theta0 cannot be drawn from its prior through any API. Liftable R-side via
setState (the state carries variance.vars/values/sizes/flags). Cost ~100-150
harness lines - build valid variance-forest structures in the flat state format,
draw leaves from the calibrated inverse-chi-square (heteroscedastic.md 3.4),
install, read the variance channel back - plus the moment self-check that the
installed s^2 matches the leaf prior, the load-bearing part. No engine change;
roughly a day.

## Exact-gate census - checking the memo's "14 of 16 single-tree"

The gates live in .github/workflows/exact-gates.yaml (push + PR on bartcore);
there is no benchmarks/exact-gates.yaml. The workflow runs 20 scripts. The
census's denominator of 16 is the 3 balance gates plus the 13 `-exact` gates.

    m = 1 (14): bd/change/swap-balance, aft-, bcf-, bcf-weak-, bcf-restricted-,
       categorical-, linear-, hazard-, hurdle-, negbin-, ordinal-, t-exact
    m > 1 (2): heteroscedastic-exact (20 mean + 1-2 variance),
       multinomial-exact (50 in arms 1/5/6, 75 in arm 7)

The count is CORRECT but generous to the two exceptions: both exceed one tree
only where the predictor is CONSTANT - heteroscedastic-exact partA is
`x <- matrix(0.5, n, 1L)  # constant: neither forest can split`, and
multinomial-exact's multi-tree arms are intercept-only (which is what makes
their quadrature closed-form). Every tree is root-only, so tree-structure MCMC
and the backfit residual bookkeeping are gated at m = 1 in 16 of 16, not 14 of
16; the extra trees exercise only the ensemble leaf/level draws.

The 4 scripts outside the 16: hazard-reduction (40 trees) and hurdle-reduction
(30) run at ensemble size but their oracle is a hand-expanded ENGINE fit,
bitwise - not an independent posterior; monotone-reference is m = 1;
logistic-reference part 1 (the gate) is m = 1 and its part 2 at ntree = 50/200
against BART::lbart/pbart is the ONLY true ensemble-scale comparison against an
independent implementation on the push path - and it does not gate (anyFailure
is fixed before part 2 runs) and is skipped unless BART is installed.

## What still has no ensemble-scale calibration evidence

After this run: aft and hetero, both with their own sampling code and both gated
only at m = 1 (or m = 20 under a constant predictor). Every other feature-matrix
row now carries an ensemble-scale SBC verdict, directly or by a bitwise
reduction to an arm that does. The pair is the honest answer to "which families
are trusted from a single-tree derivation only".
