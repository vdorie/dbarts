# Creation-time leaf-scale pin (the FIX-A shape) - scoping memo

Read-only, tip b102e17c. Item: consolidated-report.md (e) S10 (i) / (h), from
calibration-sbc.md's "not covered" pair; priced by runsbcbcf-repair.md at
~30-60 engine lines.

VERDICT: the capability already shipped. `nameable-calibration`
(docs/design/nameable-calibration.md, creation half c2a7e89b, mid-chain half
landed) IS the FIX-A pin, on the R surface, live for gaussian, aft and grouped.
The item is a harness change - zero engine lines, zero draws, zero consumer
impact - probed on all three, demonstrated on the grouped `sigma` artifact.

## 1. What "pin at creation" means, precisely

Two things are called the leaf scale; only one was ever unpinned.

- INTERNAL `forest.leaf.scale` ([[chain.hpp:633@658869ac]]) is computed ONCE at
  construction, `resolvedNodeScale(nodeScale, priorScale) / sqrt(m)`, and NO
  mutation path recomputes it - setResponse/setOffset/setData/setWeights never
  write it. Already pinned at creation, and always was.
- IN FORCE, response units: `prior.scale = leaf.scale * fitScale * sqrt(m)`
  ([[chain.hpp:3973@658869ac]]), `fitScale` the transform's multiplier - `range_` under
  gaussian/t/aft/grouped ([[model.hpp:2920@658869ac]], 3926), 1 under the latent families.
  What moves is the TRANSFORM, not the leaf scale.

So the problem was never "the engine recomputes a scale mid-run": two
INDEPENDENTLY CONSTRUCTED samplers get different response-unit leaf priors, each
deriving `range_` from its own `y`, so a harness that must rebuild per
replication (aft: censoring status is structural, [[model.hpp:3868@658869ac]]; grouped:
setData refused, [[R_interface_bartcore.cpp:4634@658869ac]]) fits under a different prior
from the one theta0 came from. `resolvedNodeScale` ([[chain.hpp:3963@658869ac]]) is the pin:
a finite `priorScale` is divided BY the transform, so the response-unit prior is
the named constant whatever `y` built it.

## 2. What each path does today (probed, not read off)

`sbc-logs/leafscale-pin-probe-calibration.R`, private lib at HEAD; anchor y
range 5, "rebuild" y range 20 (same rows, `4y + 7`).

    channel                          prior.scale   prior.mean
    creation, default node prior       2.5 / 10.0    0.0 / 7.0
    creation, normal(2, scale = 2.5)   2.5 /  2.5    0.0 / 7.0
    setResponse(updateScale = FALSE)   unchanged     unchanged
    setResponse(updateScale = TRUE)    2.5 -> 10.0   0.0 -> 7.0
    setData(new data)                  2.5 -> 10.0   0.0 -> 7.0
    setData under a NAMED scale        2.5 -> 10.0   (intent survives on the
                                                     model; setModel re-issues)

Matches nameable-calibration.md section 4. Two halves: `prior.scale` is
pinnable, `prior.mean` (the shift, `range_*0.5 + min_`) is NOT - it stays
data-derived under a rebuild, the runsbcbcf-repair survey addendum's "a correct
override must pin scale AND shift", still true. The shift's lever is the offset
channel (dbartsSampler-class.Rd `prior.mean` item), pulled AFTER construction:
the transform is taken from `y - offset`, so a creation-time offset cancels.

Family reach, probed (`leafscale-pin-probe-aft-grouped.R`): the named scale pins
across rebuilds for AFT (internal bartcore surface, status on the
`bartcore.survival` control attribute) and GROUPED - 2.5/2.5 named against
2.5/10.0 default in both - and `$setCalibration(prior.scale =)` writes mid-chain
on both, with `setOffset(rep_len(-prior.mean, n), updateScale = FALSE)` accepted
on aft. Only BCF and multinomial refuse a named scale; neither is in scope.

## 3. User-visible behavior

Under option B: NONE. Under A it depends on which channel a new pin covers (a
creation-only one duplicates `node.prior = normal(scale =)`):

- A pin held through `setData` WOULD change the model users get: `setData`
  re-anchors today (section 2), and [[dbartsSampler-class.Rd:419@658869ac]] plus dbarts.Rd's
  "Naming the leaf calibration" both PROMISE it re-anchors while leaving the
  named intent alone.
- `setResponse(updateScale = TRUE)` keeps rescaling either way - its documented
  purpose ([[dbartsSampler-class.Rd:144@658869ac]], "should only be TRUE during burn-in") -
  and is already refused where a calibration is stated against a transform
  never restated (multi-forest, heteroscedastic, grouped under gaussian/t/aft:
  `refuseGroupedScaleUpdate`, [[R_interface_bartcore.cpp:2723@658869ac]]). Untouched here.

Model-change risk from "pinning": nil under B and nil under a default-off A -
except that A's only non-duplicate content IS that setData arm.

## 4. Draw impact and the equivalence scenarios

- Option B: no shipped code, no draws.
- Option A, default-off creation pin: bitwise-neutral. No equivalence scenario
  names a `prior.scale` (benchmarks/R/equivalence.R greps clean) and
  `setForestPriorScale` already skips a write reproducing what is in force
  ([[chain.hpp:1208@658869ac]]).
- Option A extended to `setData`: DRAW-MOVING on exactly one canonical
  scenario, `setdata` ([[equivalence.R:139-146@658869ac]], the only `$setData` caller;
  gaussian, n 400 -> 500, so the transform genuinely moves) - one re-record
  plus any tinytest snapshot replaying a setData path. NO scenario calls
  `setResponse`; the other mutation scenarios drive setPredictor / setWeights /
  setActiveRows / setTestOffset, none of which re-anchor, and
  bcf-equivalence.R / multinomial-equivalence.R never mutate the response.
- Consumer (inst/include/dbarts/dbarts.h): none under B; under A, none if it
  rides fields the reshape already carries - but nameable-calibration's flat-C
  half is already an item inside dbarts.h reshape S1, so A competes with it.

## 5. Live demonstration - the grouped `sigma` artifact

The lane's driver (`sbc-logs/run-arm.R`) calls `runSbcGrouped` at the default
`swap = FALSE`, the REBUILD arm - which is why the artifact shows. Three arms,
one config (grouped gaussian, n = 160, 8 groups, rel.scale 0.2, 20% zero
weights), m = 50, R = 120, L = 150, thin = 30, burn 2400 sweeps, 5% band 0.1172
(logs and rank matrices under `sbc-logs/leafscale-pin-grouped-*` recording
not retained).

    arm                                sigma ecdf  chisq p  verdict  worst other
    rebuild, unpinned (the recorded)     0.2360    0.000    FLAG     0.0884 PASS
    rebuild, scale pinned (option B)     0.1052    0.252    PASS     0.1165 PASS
    swap arm, scale + shift pinned       0.0680    0.478    PASS     0.1048 PASS

- The pin RETIRES the artifact: 0.2360 -> 0.1052, chisq p 0.000 -> 0.252, all 10
  functionals inside the band. The change is two lines - read the anchor's
  `prior.scale` off a sampler built on `yBuild`, then `cfg$nodePrior <-
  dbartsPriors$normal(k, scale = S)`; the generator is bitwise untouched.
- The residual is the SHIFT half: in the scale-pinned arm four cells' chisq p
  fall to 0.003-0.037 (avg.f, b1, b2, f.star5) with ecdf inside the band, while
  the swap arm - pinning scale AND shift by construction (one fit,
  `setResponse(updateScale = FALSE)`) - is clean everywhere, chisq p 0.071-0.784.
  The unpinned rebuild's in-force prior.scale runs ~2.3x the anchor's; its
  unpinned prior.mean (mid-range of y0) is second-order.
- Zero-cost today: `sbc.R grouped-gaussian-swap` already exists and pins both;
  rebuild-plus-pin is what aft needs, being unable to swap.

## 6. Options

A - the pin as FIX-A describes (engine). Buys nothing not already shipped
unless it extends to `setData`. Costs ~30-60 engine lines plus an rchk note,
Opus tier, duplicate public surface, a second entry competing in reshape S1.
Draws neutral default-off, draw-moving on `setdata` if it pins that path;
consumer prose at worst. Risk: that arm contradicts two shipped Rd promises.

B - harness only (RECOMMENDED). Buys an ensemble-scale aft SBC arm (rebuild per
replication under a named scale; aft and hetero are the last families with own
sampling code and no ensemble-scale evidence) and retires the grouped sigma
artifact, demonstrated above. Costs ~2 lines for the grouped arm, ~80 for the
aft arm (calibration-sbc.md's estimate, unchanged), ~10 for the shift recipe;
Sonnet under an Opus read, benchmarks/R only. No draws, no consumer impact, no
model risk: it names a prior the shipped API already lets a user name.

C - leave as is. Buys nothing; aft ships with calibration evidence at m = 1
only and the grouped `sigma` FLAG stays as an artifact to be re-explained at
every re-run. Its aft alternative - a censoring-status setter - is NEW PUBLIC
SURFACE (~100-150 lines, three layers), landing before 1.0-0 or not at all.

## 7. Recommendation

OPTION B, pre-RC. DECIDING FACT: the pin FIX-A was priced to build is already
on the shipped surface as `node.prior = normal(scale =)` / `prior.scale`, live
on aft and grouped - probed at HEAD, and shown to convert the grouped `sigma`
FLAG (0.2360, p 0.000) to PASS (0.1052, p 0.252) through a two-line harness
edit. S10 (i) is a benchmarks item, not an engine item.

Riders: (1) pin the SHIFT as well as the scale on any rebuild arm, via
`setOffset(rep_len(anchorMean - prior.mean, n), updateScale = FALSE)` after
construction - section 5's third row is the evidence it matters; (2) the aft arm
needs nothing else from the engine - it draws sigma conjugately like gaussian
([[chain.hpp:648@658869ac]]), its reported-scale sigma prior is range-independent under an
explicit `sigma =` ([[sbc.R:38-42@658869ac]]), and a zero-censoring control reduces bitwise
to the gaussian arm ([[survival.md:53-56@658869ac]]), the self-check to run first.
