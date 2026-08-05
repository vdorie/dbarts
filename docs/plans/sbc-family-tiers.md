# sbc-family-tiers

agent: opus implementer; fable adjudicates the rank histograms
rng: neutral - benchmarks/ and .github/ only; gaussian CI ranks replay exactly
budget: one sbc.R block per family, a discrete-rank helper, a CI matrix

## Goal

A recorded SBC verdict for every shipped family, or the reason it is ill-posed.

## Context

- sbc-calibration.md; sbc-ci-gate.md, whose admission rule step 4 replaces; the
  design notes. runSbcBCF errors at HEAD (9030d93): OUT of scope, fixed separately.

## Decision - scope

The fork: which families, at which prior. IN below at R=200, L=150, thin=30, burn
MEASURED (step 2). Changes it: a k=8 cost over budget makes nbinom manual-only; a
burn ladder against the 72000 floor re-opens burn for every arm.

    R=200 at the SHIPPED 72000-sweep burn, 38/48/56 us/sweep: 10-15m, not ~2m
    ordinal   gamma_2..K-1, eta, agg p  marginal-MH cutpoints     ~14m
    nbinom    r, avg mu, agg psi [k=8]  grid r cond., PG y+r         ?
    t         sigma, nu, avg f, agg f*  lambda mixture, nu grid   ~12m
    multinom  agg p_k(x*), raw f_ik     interleaved PG, centering    ?

- nbinom at a TIGHTENED k = 8 (psi sd = node.scale/k = pi*sqrt(3)/8 ~ 0.68, vs 2):
  simulatePolyaGammaShape loops sum(y_i + r) times/sweep, lognormal-tailed and
  unbounded under default-k draws (13.5x over six, tail reps hours), so it is not
  budgetable; a tightened prior still validates NB. RE-MEASURE the sweep cost at 8.
- multinom also ranks raw f_ik at a few (i, k) cells via bartcoreForestFits (the
  BCF accessor): the softmax flat direction is no prior symmetry, so f_ik SBC is
  well-posed and tests MultinomialForestCombiner::afterCombine's level-centering
  draw, whose precision (sum of invV*n) approximates the exact per-leaf
  conditional - the pre-registered suspect if f_ik alone flags.
- r and nu are DISCRETE, state-only: rank = #{draws < theta0} + Uniform{0..#ties}
  - theta0 takes its own iid tie-breaker beside the L draws, so the null is
  rankUniformity's uniform on {0..L}; collect by per-sample run(0,1)+storeState.

OUT, with the change that would enable each:
- aft: status is fixed at creation, so a status-varying DGP forces a rebuild that
  re-anchors range(y0). Enabler: a status setter.
- hazard: the at-risk person-period expansion makes the design - and so the tree
  prior (sampleTreesFromPrior's cut grid, empty-node collapse) - depend on y0,
  breaking exchangeability; second, it owns no sampling code (hazard-reduction.R
  gates it bitwise vs a hand-expanded fit).
- hurdle: no new sampling code; same design-depends-on-y0 defect (y > 0 subset).
- heteroscedastic: prior draws never reach varianceForest_; liftable R-side via
  setState (state carries variance.vars/values/sizes/flags). Deferred, not blocked.
- monotone: the enabler (a constrained prior draw) LANDED 173a710 the day this
  plan was written (monotone-prior-draw.md), so monotone is now liftable; first
  follow-on arm after the initial matrix, not in it.

## Constraints

R entry points of the installed package only, no engine edits. Out of scope: BCF
repair, zero-inflation, grouped families, xbart, Tier A-C, sampled-k/DART/tau,
and the critique's monotone bonus finding (a SEPARATE item, not this plan).

## Steps

1. Discrete tie-break helper + a sbcCheckSigmaPrior-style self-check: a synthetic
   discrete conjugate case with closed-form posterior must pass rankUniformity.
2. Measure a per-family burn ladder first - 72000 was a BCF-specific floor - and
   for ordinal pre-commit to the A4e chain-length ladder before reading any
   cutpoint/eta flag as a defect (ordinal.md section 9's cutpoint-shift ridge).
3. Per family in table order: prior draw + moment check, likelihood, functionals,
   thin/burn, R=200 run, verdict. ordinal and nbinom REBUILD per rep (safe: scales
   are 1, node.scale constant; it clears the kept gamma_/r_ that would break rank
   iid-ness); only t reuses a pinned sampler (setResponse cold-inits nu, lambda).
4. sbc.yaml matrix (fail-fast: false, timeout from step 2). Admission: the alpha
   Bonferroni'd to 0.05/M over the matrix's functional count (full-matrix pass
   ~0.95 on a fresh stream); drop gamma_1 (pinned at 0), rank gamma_2..K-1 with
   K >= 4, aggregate the (k, x*) cells. Seed pinned; risk: a draw-shifting commit
   or CI image re-rolls the stream, so the gate stays non-blocking.

## Verification

- `Rscript benchmarks/R/sbc.R gaussian 100 200 30` -> identical ranks.
- `... ordinal|nbinom|t|multinom 200 150 30`, one each -> every functional PASS.
- `... discrete-selfcheck` -> rankUniformity PASS on the conjugate case.
