# variance-forest-mixture-reinstall: setControl's mixture and perturb width reach the variance forest

Status: PLANNED 2026-10-11.

agent: opus implementer, one; opus reviewer who runs the mutants below.
rng: NEUTRAL, bit for bit, for every sampler that never calls setControl with a changed proposal.probs or
n.perturb.cuts, and for every single-forest sampler whatever it calls. SHIFTING for a heteroscedastic
sampler from such a call on: a posterior-preserving kernel change, the draws move and the stationary
distribution does not. No recorded scenario or snapshot file makes that call.
window: before 1.0-0; engine slices stay serial.
budget: ~250 lines (engine ~15, tests/cpp ~90, tinytest ~90, help 2, records ~60).

## Goal

After `setControl` with a changed `proposal.probs` or `n.perturb.cuts`, both forests of a heteroscedastic
sampler propose under the new mixture and width, as both do when the sampler is created under that control.
A live sampler and its own copy or reload run the same kernel. The tier is "Changes draws"
([Process by risk](README.md#process-by-risk)).

## Context

Facts, read and run on d5b637ee (shipped build, 2026-10-11):

- At creation the variance forest takes the control's mixture and width:
  [`buildVarianceForest`](../../src/bartcore/chain.hpp) copies the seven fields of
  [`SamplerOptions`](../../src/bartcore/chain.hpp) the mean forest is given, which the bridge's
  [`parseProposalProbs`](../../src/R_interface_bartcore.cpp) reads off the control.
  [`varianceForest`](../../R/model.R) takes `vars`, `n.trees`, `base` and `power` only, no mixture or
  width of its own, and [`sweepVarianceForest`](../../src/bartcore/chain.hpp) states the intent: "the
  variance forest reads the same mixture as the mean forest, so an all-zero one freezes both". Run: under
  a mixture frozen at creation the variance surface holds one cell at every one of 30 sweeps.
- [`dbartsSampler$setControl`](../../R/dbarts.R) stores the control and, where the mixture or the width
  moved, calls `$setModel(model)`; the bridge re-reads the control's mixture and
  [`Chain::setModel`](../../src/bartcore/chain.hpp) writes it to the first forest alone. Run: frozen by
  setControl after 50 growing sweeps, the mean trees stand at all 30 following sweeps and the variance
  surface's cells move at all 30; setControl before any run to a perturbing mixture at width 3 draws
  bitwise what creation under it draws on a mean-only sampler and not on a heteroscedastic one.
- A copy or a reload re-creates the engine from the stored control, so its variance forest takes the new
  mixture. Run: after that frozen setControl, the live sampler's variance cells move at 30 of 30 sweeps
  and its copy's at 0 of 30. The live sampler and its copy run different kernels today.
- The tree-move probabilities are not model and stay with setControl (dec-B254). The help's two sentences
  that the variance forest keeps what it was created with record the defect, from a review, not a ruling.
- Other forests. A sampler whose forests carry amplitudes (BCF-shaped) and a multinomial one give every
  forest the control's mixture and width at creation
  ([`ForestStructureSpec`](../../src/bartcore/combiner.hpp),
  [`MultinomialForestSpec`](../../src/bartcore/combiner.hpp)), and both refuse the install whole:
  `$setModel` stops in [`refuseAmplitudeMutation`](../../R/bartcore.R) or
  [`refuseCountsMutation`](../../R/bartcore.R), the bridge behind them in
  [`refuseMultiForestMutation`](../../src/R_interface_bartcore.cpp), and setControl rolls the stored
  control back (run: both do). A hurdle fit is two single-forest samplers, each with its own control, and
  takes no variance forest. The variance forest is the only second forest on a sampler that accepts the
  install.

## The rule

`setControl(control)` means the sampler now runs under that control. Every forest that took a setting from
the control when the sampler was created takes the new value when the control changes; a forest whose
setting is its own keeps it. The mixture and the perturb width are the control's on every forest; the
variance forest's `base` and `power` are its own and are not touched here.

## Constraints

- No signature change in `inst/include/dbarts/dbarts.h` (the flat C API installs no mixture), the facade or
  the bridge; no R code change. Multi-forest samplers keep refusing the install.
- No NEWS: the variance forest is new in 1.0-0, so no released version reaches the defect.
- Out of scope: the variance forest's `base` and `power` at `$setModel` (dec-B254: still to be built).

## Steps

1. [`Chain::setModel`](../../src/bartcore/chain.hpp): where a variance forest exists, write the six
   probabilities and the width to it beside the mean forest's. Nothing else on the variance forest moves:
   its tree prior, trees and factors stand, and the scale leaf is recalibrated as today.
2. Help: the `proposal.probs` and `n.perturb.cuts` items of
   [man/dbartsControl.Rd](../../man/dbartsControl.Rd) say the install reaches both forests of a
   heteroscedastic sampler. TODO item removed.

## Tests

- tests/cpp, a new test in test_sampler.cpp, reading the installed fields through a new
  [`TestPeer`](../../tests/cpp/test_peer.hpp) accessor: after `setModel` the variance forest's six
  probabilities and width are the model's and its base and power are creation's; a heteroscedastic sampler
  handed the model before any run draws bitwise what creation under it draws, as a mean-only one does.
- tinytest, a new file test-variance-forest-mixture.R, heteroscedastic single chain: frozen by setControl
  after growing sweeps, the mean trees (getTrees) and the variance surface's cells (rows grouped by equal
  getVariance values) both stand over the following sweeps while the surface's values keep moving;
  setControl before any run equals creation under the control, once for a changed mixture and once for a
  changed width alone; a live sampler and its copy both stand after the frozen install; a sampler handed
  its own unchanged control draws what a twin that never calls setControl draws.
- Mutants, each failing a named expectation: the block removed; the width left out; one probability left out.

## Verification

- `R CMD INSTALL --preclean` into the slice's libraries (plain and reference), `buildInfo()$mode` checked;
  tests/cpp plain and under address and undefined sanitizers; the full tinytest suite, one process,
  at_home, 0 failures.
- Equivalence trio on the reference build, `compare --bitwise --strict-coverage`: equivalence-308c45be 55
  of 55 identical, bcf-equivalence-1b7d730c 15 of 15, multinomial-equivalence-80b1c8d4 11 of 11, 0
  skipped. A scenario that moves stops the slice.
- The four seeded snapshot files pass unchanged on the reference build.
- Exact gates in quick mode that fit a variance forest: heteroscedastic-exact.R, aft-exact.R, aft-hetero-pit.R.
- `lintr::lint_package()`, `air format --check .`, and the three tools/ checks, each on its own exit
  status; `R CMD check --as-cran` from a tarball built from a clean copy staged outside the tree.

## Calls made

- The fix is the reinstall, not a record of why not: creation gives the variance forest the control's
  mixture and width, and a copy already takes the new ones.
- The write lives in `Chain::setModel`, beside the mean forest's, so `$setModel` alone re-installs the
  stored control's mixture on both forests: the values already in force, a bitwise no-op.
- `Chain::setModel` is not widened to every mean forest: the R surface and the bridge refuse the install
  on a multi-forest sampler, and the rest of the install (leaf scale, tree prior) is single-forest there.

## Not verified

- No exact-posterior gate changes the mixture part way through a run on a heteroscedastic sampler. The
  claim that the stationary distribution stands rests on each mixture being a valid kernel for the
  variance forest, which creation under any mixture already relies on.
- Several chains: by construction ([`Sampler::setModel`](../../src/bartcore/sampler.hpp) loops them), no test.
