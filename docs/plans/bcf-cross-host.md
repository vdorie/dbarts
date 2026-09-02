# BCF equivalence baseline: cross-host mode (D7)

Status: LANDED 2026-08-26 at 3f532af2 (design record cbb1cb97; see the
landing note). Read-only; no repo edits. Revised after an independent
blind critique and the coordinator's adjudications - sections 1, 2, 4, 5,
7 and 8 carry the rulings.

Spec: docs/plans/prerc-surface-freeze.md D7 ([[docs/plans/prerc-surface-freeze.md:72-80@cbb1cb97]]) and its Sequencing
line ([[docs/plans/prerc-surface-freeze.md:100@cbb1cb97]]). TODO `bcf-baseline-cross-host` ([[TODO:33-34@cbb1cb97]]), which absorbs
the cross-host half of `equivalence-harness-statistical-mode`
([[TODO:66-73@cbb1cb97]]). D7 is the last open pre-RC item (`rc-gate`, [[TODO:164-165@cbb1cb97]]).

Three facts settled by reading, each of which changes the shape of the
work:

- D7's seven-name exemption list is WRONG on two names (section 1). The
  correct list is SIX, and one of D7's names (`varcount`) is a DRAWS-axis
  channel - exempting it would gut the gate.
- The current `bcf-equivalence-6e3b9fb8.rds` DOES carry enough to gate
  cross-host (section 3): it stores every raw draws-axis channel, and the
  absent `summaries` field is a pure function of what is stored. The
  RC-tip re-record is a refresh, and the box validation runs this week
  against the file in the tree.
- A single |z| bar cannot be the cross-host gate (section 2). bcf's z is
  over 40 AUTOCORRELATED draws of ONE chain, not over independent seeds
  as equivalence.R's is, and at the measured autocorrelation the |z| = 4
  bar tolerates a shift of ~0.9 posterior sd per cell. The verdict is
  TWO-TIER: a tight deviation bound that is the real gate, and a
  decoupled, explicitly-labelled statistical tier beneath it.

Budget: ~90 lines in benchmarks/R/bcf-equivalence.R, ~85 in
benchmarks/R/multinomial-equivalence.R, ~12 in benchmarks/README.md, ~10
in .github/workflows/exact-gates.yaml, 2 tokens in equivalence.yaml. No
R/, src/, inst/ or man/ file moves; no baseline re-records at this slice.

## 1. Channel taxonomy, verified against the script and the baseline

`recordChannels` ([[benchmarks/R/bcf-equivalence.R:412-424@23b9cde7]]) writes seven
channels plus a derived `summaries`; five scenarios append a verdict at
the call sites ([[benchmarks/R/bcf-equivalence.R:291@23b9cde7]], [[benchmarks/R/bcf-equivalence.R:312@23b9cde7]], [[benchmarks/R/bcf-equivalence.R:338@23b9cde7]], [[benchmarks/R/bcf-equivalence.R:369@23b9cde7]], [[benchmarks/R/bcf-equivalence.R:392@23b9cde7]]). Nine names in all. Shapes
and storage modes below are the recorded file's, not inferred:

| channel | shape | type | source | draws axis | kind |
| --- | --- | --- | --- | --- | --- |
| `mu` | 200 x 1 | double | `bartcoreForestFits(bc, 0)` [[benchmarks/R/bcf-equivalence.R:414@23b9cde7]] | no | continuous |
| `tau` | 200 x 1 | double | `bartcoreForestFits(bc, 1)` [[benchmarks/R/bcf-equivalence.R:415@23b9cde7]] | no | continuous |
| `glue` | 3 x 1 | double | `bartcoreForestAmplitudes` [[benchmarks/R/bcf-equivalence.R:416@23b9cde7]] | no | continuous |
| `sigma` | 40 | double | `result$sigma` [[benchmarks/R/bcf-equivalence.R:417@23b9cde7]] | YES | continuous |
| `train` | 200 x 40 | double | `result$train` [[benchmarks/R/bcf-equivalence.R:418@23b9cde7]] | YES | continuous |
| `varcount` | 4 x 2 x 40 | integer | `result$varcount` [[benchmarks/R/bcf-equivalence.R:419@23b9cde7]] | YES | combinatorial |
| `varcount.tau` | 4 x 1 | integer | `bartcoreForestVariableCounts(bc, 1)` [[benchmarks/R/bcf-equivalence.R:420@23b9cde7]] | no | combinatorial |
| `accepted` | 1 | logical | 3 scenarios | no | combinatorial |
| `installed` | 200 | logical | 2 scenarios | no | combinatorial |

TWO axes, and both matter. The DRAWS-axis column decides what is
exempted (D7's decision). The KIND column decides how a gated channel is
compared cross-host, and it is the split the single-|z| design missed: a
Welch z on an integer channel is a bitwise test wearing a z. `varcount`
in `restricted` has 2 structurally zero-variance cells (the moderator
restriction to {x1, x3} leaves x2 and x4 unsplit in the tau slab), which
give 0/0 -> NaN when the two runs agree and diff/0 -> Inf when they do
not; the other 11 scenarios have none. Measured across the recorded file:
the printed summary count is 209 per scenario (1 sigma + 200 train + 8
varcount), 2508 over the 12, uniform - it is the count of INFORMATIVE
cells that is not.

Draws-axis (gated) is exactly `statChannels` ([[benchmarks/R/bcf-equivalence.R:113@23b9cde7]]): `sigma`, `train`,
`varcount`. Snapshot (exempt cross-host) is the other six: one
point-in-time query of the live sampler after the last sweep.

Why the continuous snapshots can never pass bitwise off-host: `makeData`
([[benchmarks/R/bcf-equivalence.R:89-97@23b9cde7]]) builds `y` through `sin(pi * x[, 1])`, and `sin` is platform
libm, so the INPUT DATA differs in the last ULP macOS-vs-Linux before any
engine code runs.

Against D7's list `(mu, tau, glue, varcount, forestFits, accepted,
installed)`:

- `forestFits` is NOT a channel of this harness. It is
  multinomial-equivalence.R's name for the same quantity ([[benchmarks/R/bcf-equivalence.R:119-122@23b9cde7]]); bcf
  records the two forests separately as `mu` and `tau` ([[benchmarks/R/bcf-equivalence.R:130-131@23b9cde7]]). Drop
  it - it duplicates two names already on the list, imported from the
  sibling script.
- `varcount` is DRAWS-axis, 4 x 2 x 40, in `statChannels` since the
  forest-axis widening (a16d703c), and the only per-draw guard on BOTH
  forests' tree structure ([[benchmarks/R/bcf-equivalence.R:20-24@a16d703c]]). The intended name is `varcount.tau`,
  the 4 x 1 cumulative final-state query at [[benchmarks/R/bcf-equivalence.R:136@a16d703c]] - which the TODO's own
  garbled parenthesis ([[TODO:84@a16d703c]] pre-collapse) listed on both sides of its
  semicolon; gone now that TODO's equivalence-harness-statistical-mode
  entry is collapsed to its closed summary.

Corrected exemption list, six names:

    snapshotChannels <- c("mu", "tau", "glue", "varcount.tau",
                          "accepted", "installed")

Make the D7 error unrepresentable. Two asserts, and the second must sit
in the COMPARE LOOP over `names(a)`, not in `recordChannels`: `accepted`
and `installed` are appended at the call sites, outside `recordChannels`,
so an assert placed there cannot see them and a new verdict channel added
the same way would slip through.

    # at script load
    stopifnot(length(intersect(snapshotChannels, statChannels)) == 0L)
    # inside the compare loop, per scenario
    unclassified <- setdiff(names(a),
                            c("summaries", statChannels, snapshotChannels))
    if (length(unclassified) > 0L)
      stop("unclassified channel(s) in ", name, ": ",
           paste(unclassified, collapse = ", "))

The second is the one that matters going forward: a channel added later
falls into neither list and STOPS the compare, forcing a taxonomy
decision instead of silently defaulting to gated-and-red off-host.

## 2. The flag, and the two-tier cross-host verdict

Spelling `--cross-host`, matching `--strict-coverage`
([[benchmarks/R/equivalence.R:34-35@a16d703c]]). Parsed position-independently and
removed from `args` before the mode is read, so it composes with `quick`
and cannot be mistaken for a mode or a path:

    crossHost <- "--cross-host" %in% args
    args <- setdiff(args, "--cross-host")

Compare-mode only; `record` is host-local by definition and says so if
passed the flag. Landed as `compareCrossHost` ([[benchmarks/R/equivalence.R:224-226@a16d703c]]), called from the
compare loop only when `crossHost` is set ([[benchmarks/R/equivalence.R:793@a16d703c]]):

    exempt <- intersect(snapshotChannels, names(a))
    gated <- setdiff(names(a), c("summaries", exempt))

### 2.1 Why one |z| bar cannot be the gate

equivalence.R's z is over 20 INDEPENDENT SEEDS ([[equivalence.R:1433@a16d703c]]:
`colMeans` of a seeds x summaries matrix). bcf's `drawSummary`
([[equivalence.R:118-123@a16d703c]]) reduces over 40 AUTOCORRELATED draws of a single chain, and
`n` in the Welch denominator ([[equivalence.R:514@a16d703c]]) is that 40. Measured on the recorded
`default` scenario's 200 train cells: lag-1 autocorrelation median 0.423,
AR(1) effective sample size median 16.2 (min 3.4), so the denominator
understates the true MCSE by sqrt(40 / 16.2) = 1.57x.

Vacuity, quantified. At the nominal n = 40 the |z| = 4 bar tolerates a
per-cell mean shift of 4 * sqrt(2/40) = 0.894 posterior sd: 0.108
absolute on `train` (median posterior sd 0.121) and 0.0086 on `sigma`
(sd 0.00965). With the ESS correction the bar is looser still -
4 * sqrt(2/16.2) = 1.41 posterior sd, 0.170 on `train`, which is 0.85 of
the residual sd (sigma ~ 0.1996). A defect that shifts every fitted value
by four fifths of the residual scale passes. The statistic is also
bimodal, so the "1-4 band" is a band it essentially never occupies:
stream locked -> |z| ~ 1e-13; one flipped MH decision -> decoupled chains
where the bar is a true ~2.5 over 209 cells per scenario, red with near
certainty. The bar is either vacuous or a coin flip; it is not a gate.

### 2.2 Tier 1 - the gate (LOCKED)

Cross-host under a locked RNG stream the two builds agree to rounding.
That is what the gate tests, per channel per scenario:

- CONTINUOUS gated channels (`sigma`, `train`; any continuous draws-axis
  channel added later): every cell within `abs(a - b) <= atol + rtol *
  abs(a)`, `rtol = 1e-8`, `atol = rtol * max(abs(a))` over that channel
  in that scenario. Report the worst ratio
  `max(abs(a - b) / (atol + rtol * abs(a)))`; PASS iff <= 1. On `train`
  that atol is 5.7e-8 and on `sigma` 2.2e-9 - about 1.8e6 times tighter
  than tier 2's 0.108, and non-vacuous against every defect class the
  harness exists to catch.
- COMBINATORIAL gated channels (`varcount`): `identical()`. Split counts
  are integers; under a locked stream they reproduce exactly, and the
  0/0 -> NaN and diff/0 -> Inf cells of section 1 stop being a bitwise
  sub-gate hidden inside a statistical one.

CALIBRATION CLAUSE: `rtol = 1e-8` is a proposal, not a measurement. The
prior cross-host evidence (`~/equiv-3a977b6d.log`, 37 scenarios, 2183
summaries, every line `max |z| = 0.00`, exit 0) bounds the MEAN shift,
not the per-draw deviation, so it cannot fix the constant. The box run
(section 4) REPORTS the observed worst ratio; if it exceeds 0.01 - within
100x of the bound - widen `rtol` to 100x the observed maximum and record
the measured value here and in the MANIFEST row. Do not ship a bound the
first real run sits just under.

### 2.3 Tier 2 - decoupled, and weak by construction

Reached only when tier 1 fails. Its job is adjudication - did the
posterior move, or did the stream merely decouple - never certification.

- Continuous channels: Welch z with an ESS-AWARE denominator. Per cell,
  `ESS = n * (1 - r1) / (1 + r1)` from the lag-1 autocorrelation of that
  cell's draws (pooled over the two runs, floored at 2), and
  `z = (mean_a - mean_b) / sqrt(var_a/ESS_a + var_b/ESS_b)`. This is
  CONSERVATIVE relative to the nominal statistic: it enlarges the
  denominator by ~1.6x at the measured autocorrelation, so a mere
  decoupling is less likely to be called a posterior move. Bars stay
  |z| > 3 warn, |z| > 4 fail.
- Combinatorial channels: reported DIAGNOSTIC ONLY (cells differing,
  worst absolute difference). No z, no verdict - a Welch z on integer
  counts is a bitwise test in disguise and tier 1 already ran it.

The doc, the log line and the README all say the same thing about tier 2:
at the measured ESS the |z| = 4 bar tolerates a per-cell shift of 1.41
posterior sd - 0.170 on `train`, 0.85 residual sd. Tier 2 passing is not
evidence the builds agree; it is evidence the failure is not gross.

A real seeds axis is the honest fix and is a POST-RC DOOR, not this
slice: it means giving each scenario `n.seeds` independent engine seeds
and recording a seeds x summaries matrix, which changes the scenario
structure, changes `settingsList()` ([[equivalence.R:424-434@23b9cde7]]), invalidates every
existing bcf and multinomial baseline through the settings guard
([[equivalence.R:458-472@23b9cde7]]), and multiplies runtime by the seed count. Additive by rule,
so it does not gate the release.

### 2.4 Non-finite guard (was a hole)

`--cross-host` would otherwise make the harness BLIND to NaN/Inf: a
non-finite `train` makes every z NaN, `max(abs(z), na.rm = TRUE)` returns
-Inf and `n.fail` is 0 ([[equivalence.R:519-528@23b9cde7]]), so the run prints `max |z| = -Inf` and
exits 0. Today `mu`/`tau` catch that through `pointMismatch` ([[equivalence.R:548-555@23b9cde7]]);
exempting them deletes the only guard.

Count non-finite values on EVERY gated channel, in both `a` and `b`,
before the tier evaluation. Any non-finite is a FAIL with its own line:

    default        NON-FINITE: train 37 cells (b), sigma 0, varcount 0 <- FAIL

Placement, so section 6 holds: run it unconditionally under
`--cross-host`, and in the default mode only inside the mismatch branch.
`identical()` treats NaN as equal to NaN, so a bitwise-clean same-host
run would otherwise gain a line it does not have today.

### 2.5 Output

Per scenario, under `--cross-host`, three lines, each prefixed by the
scenario name so the gate-counting rule can count them (memory
`orchestrator-run-mechanics` rule 9: the terminal line is never evidence;
count the per-scenario lines and require the scenario count):

    default        exempt (cross-host): mu, tau, glue, varcount.tau - 4 snapshot channels skipped by design [3 differ, max rel dev 4.4e-16]
    default        NON-FINITE: none
    default        tier 1 PASS: max dev ratio 3.1e-03 (train 3.1e-03, sigma 8.0e-04), varcount identical

and on a tier-1 failure, that line becomes `tier 1 FAIL` and is followed
by the labelled tier-2 line, which can never be mistaken for a pass:

    default        tier 1 FAIL: train ratio 4.2e+06 (cell 118), varcount 6 of 8 cells differ
    default        decoupled: statistical - 201 summaries, ESS-adjusted max |z| = 1.83 (weak bar: tolerates 1.41 posterior sd); varcount diagnostic-only, 6 cells differ, worst 3

The exempt line names every exempted channel present in that scenario (4
on the seven plain scenarios, 5 on the three `accepted` and the two
`installed` ones) and carries a NON-GATING diagnostic - count differing,
max relative deviation for the doubles, count differing for the
integer/logical ones - so a snapshot channel that moved by 1e-3 rather
than 1e-16 is visible in the log even though it does not gate.

Cross-host PASS evidence: 12 `tier 1 PASS` lines, 12 `exempt
(cross-host):` lines, zero `FAIL`, zero `decoupled:`, exit 0. The
presence of ANY `decoupled:` line is a finding to report even when the
run exits 0.

The bitwise-clean branch at [[equivalence.R:492-496@23b9cde7]] prints `identical (all %d channels:
%s)`; under the flag it must print `identical (all %d GATED channels:
%s)` so the two modes' pass lines are not interchangeable.

## 3. Baseline format: what "re-record with summaries" means

Concretely: the recorded per-scenario list carries `ch$summaries` =
`lapply(ch[statChannels], drawSummary)` ([[equivalence.R:138@23b9cde7]]) - for each of `sigma`,
`train`, `varcount` a `list(mean=, var=, n=40)` over 1, 200 and 8 cells.
`recordChannels` has written it since 3ee6f23e, so a re-record at the RC
tip produces it with no format change: the RC-tip step is `record`, not a
schema edit.

Does the CURRENT file suffice? Measured on
`benchmarks/baselines/bcf-equivalence-6e3b9fb8.rds`: `summaries` is
ABSENT from all 12 scenarios (it predates 3ee6f23e), but every RAW
channel is present at full shape, including all three `statChannels` -
`sigma[40]`, `train[200x40]`, `varcount[4x2x40]`. Both tiers read the raw
arrays: tier 1 needs nothing but them, and `drawSummary` is a pure,
deterministic reduction of exactly those arrays, so tier 2's summaries
are derivable element for element. The loud degrade at [[equivalence.R:503-510@3ee6f23e]]
("recorded before summaries") is a false negative for this file.

    aSummaries <- a[["summaries"]]
    if (is.null(aSummaries) && all(statChannels %in% names(a)))
      aSummaries <- lapply(a[statChannels], drawSummary)

Not flag-gated: strictly better in both modes, and it cannot move the
standing same-host gate, which is bitwise 12/12 and returns at [[equivalence.R:490-497@3ee6f23e]]
without reaching this branch (section 6).

So the RC-tip re-record is a REFRESH - it pins the RC engine tip and
retires the derive path for this file - not a prerequisite. Interim
validation runs against `bcf-equivalence-6e3b9fb8.rds` as it sits in the
tree.

## 4. Validation on the x86 box, this week

Box facts, verified: Linux x86_64, Ubuntu 24.04, AMD Ryzen 7 3700X, R
4.3.3 at /usr/bin/R, gcc 13.3.0, loadavg ~0.2, `~/libshim/libtirpc.so`
present. 7G RAM (~5G available) and 144G free DISK. It retires within
days; every step below is same-week. There is no `~/.Renviron` on the
box - `R_LIBS` on every call is the private-lib pattern (one lib per
role, nothing clobbered), not an override workaround.

    git archive faf1d167 | ssh dbarts-bench 'mkdir -p ~/dbarts-val-d7 && tar -x -C ~/dbarts-val-d7'

That is the whole ship step: all three baselines are tracked
(`git ls-files benchmarks/baselines/`) with no `.gitattributes`, so the
archive carries them - no separate `scp`. `.claude` is untracked (0
files), so the archive excludes `.git` by construction and `.claude` by
accident. rsync `--exclude=.git --exclude=.claude` is the equivalent when
the working tree is dirty and that is what you want shipped.

    ssh dbarts-bench 'cd ~/dbarts-val-d7 && R CMD INSTALL --preclean -l ~/rlib-d7 . 2>&1 | tail -5'
    ssh dbarts-bench 'cd ~/dbarts-val-d7 && R_LIBS=~/rlib-d7 Rscript benchmarks/R/bcf-equivalence.R compare benchmarks/baselines/bcf-equivalence-6e3b9fb8.rds --cross-host'
    ssh dbarts-bench 'cd ~/dbarts-val-d7 && R_LIBS=~/rlib-d7 Rscript benchmarks/R/multinomial-equivalence.R compare benchmarks/baselines/multinomial-equivalence-4d9a3337.rds --cross-host'
    ssh dbarts-bench 'cd ~/dbarts-val-d7 && R_LIBS=~/rlib-d7 Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-736bfb05.rds --strict-coverage'

Plus the standing x86 confirmation legs (this slice touches no C++ and no
`R/`, so they confirm the shipped tree built, they do not gate the
change):

    ssh dbarts-bench 'cd ~/dbarts-val-d7/tests/cpp && R_LDFLAGS="$(R CMD config --ldflags) -L$HOME/libshim" make && ./test_bartcore'
    ssh dbarts-bench 'cd ~/dbarts-val-d7 && R_LIBS=~/rlib-d7 Rscript -e "tinytest::test_package(\"dbarts\")" 2>&1 | tail -20'

EXPECTED VERDICTS. Cross-host is NEVER bitwise, so no scenario may report
`identical draws` or `identical (all`. Expect 12 `tier 1 PASS` lines with
worst ratios well under 1, 12 `exempt (cross-host):` lines with max rel
dev in the 1e-16 range, zero `decoupled:` lines, exit 0; multinomial the
same over 11 scenarios; gaussian 43 compared / 0 skipped, every line
`max |z| = 0.00`. REPORT the observed tier-1 worst ratio per channel -
that number is the calibration input of section 2.2 and must appear in
the landing note.

The box run is the ENTIRE interim cross-host validation. Nothing in CI
covers this until a workflow carrying the flag is pushed to bartcore
(section 5), so a probe not run here is a probe not run.

DISCRIMINATION. Four probes, all on the box in the same session.

A correction that kills two candidate probes before they are written:
`[[combiner.hpp:334-335@736bfb05]]` (`bPriorVariance`, `sdModerate`) are DEAD on the R
path. `expandForestSpecs` returns early when `spec.forests` is non-empty
([[combiner.hpp:362@736bfb05]]), and the bridge fills it unconditionally
([[src/R_interface_bartcore.cpp:2157@736bfb05]]), reading the live values out of the
R-side vector instead: `forest.nodeScaleFactor = params[3]` ([[src/R_interface_bartcore.cpp:2174@736bfb05]]) and
`amplitudePriorVariance = params[5]` ([[src/R_interface_bartcore.cpp:2175@736bfb05]]), fed by
`sd.moderate`/`b.prior.variance` at [[R/bartcore.R:787-788@736bfb05]]. Mutating the
header defaults changes nothing any scenario runs.

- D-1', ZERO-RNG, on a gated channel's own write path.
  `[[src/bartcore/chain.hpp:5290-5291@736bfb05]]`, inside `storeSample`'s
  single-location training write, replace `fitsWithoutOffset(out);` with
  the prognostic forest's own fits:

      for (size_t i = 0; i < n; ++i)
        out[i] = scale * forests_[0].totalFits[i] + shift;

  (leaving the offset add at [[src/bartcore/chain.hpp:5294@736bfb05]] in place; `totalFits` is the same
  array `results.forestFits` copies at [[src/bartcore/chain.hpp:5351-5355@736bfb05]]). Consumes no RNG,
  mutates no state, and displaces ONLY the reported `train` channel by
  the treatment contribution `b_z * tau(x)`, ~0.3-1.0 absolute. It leaves
  `mu`, `tau`, `glue`, `sigma`, `varcount` untouched - so it can only be
  caught through a gated channel, which is exactly the property under
  test. MUST fail tier 1 AND tier 2 in all 12 scenarios. Build into a
  separate lib (`-l ~/rlib-d7-mut`); after reverting, `touch` the header
  (an `mv` back preserves the mtime and the next install keeps the
  mutated object).
- D-2, FLAG IS LOAD-BEARING. The same mutated build, cross-host, WITHOUT
  `--cross-host`. Must exit 1 with `also MISMATCH in mu, tau, glue (no
  statistical fallback) <- FAIL` ([[src/bartcore/chain.hpp:551-555@736bfb05]]) on every scenario - the ULP
  divergence is unconditional, so the unflagged mode still refuses
  cross-host and the flag is not a no-op.
- D-3, TAXONOMY. Add `"varcount"` to `snapshotChannels`; the load-time
  `stopifnot` must abort. That is the D7 list error, executed.
- D-4, STREAM-PERTURBING. `[[src/R_interface_bartcore.cpp:2174@736bfb05]]`, scoped to
  the treatment forest: `forest.nodeScaleFactor = (f == 1 ? 1.2 : 1.0) *
  params[3];` (equivalently `sd.moderate = 1.2` at [[R/bartcore.R:748@736bfb05]]).
  Widens tau's node prior 20%, so the leaf draws change magnitude, the
  stream decouples, and the posterior genuinely moves. MUST fail tier 1.
  Whether tier 2 then catches it is the calibration evidence for how weak
  tier 2 is - report the ESS-adjusted max |z| per scenario either way,
  and put it in the landing note. If tier 2 passes a 20% prior change,
  say so in the doc; that is the finding, not a defect in the probe.

Cost: two `--preclean` installs (~5 min each on the Ryzen); every compare
is seconds. Report the measured wall time of each harness's compare - it
is the input to section 5's CI cost.

## 5. Durability: which CI host, and what is dormant

The correction that decides this: of the ten workflows, only
`check-standard.yaml` is on main. `equivalence.yaml` - which already
carries `bcf-equivalence` ([[R/bartcore.R:64-89@736bfb05]]) and `multinomial-equivalence`
([[R/bartcore.R:91-116@736bfb05]]) jobs on `ubuntu-latest`, i.e. x86_64 Linux, against
arm64-recorded baselines - is ABSENT from main, and `workflow_dispatch`
exists only for default-branch workflows. It is therefore neither
scheduled nor dispatchable today. "Dispatch it by hand" is not available;
the cross-host CI half is DORMANT until bartcore merges.

`exact-gates.yaml` is push-triggered on `[main, master, bartcore]`
([[R/bartcore.R:28-35@736bfb05]]), does not ignore `benchmarks/**`, and a push run executes the
PUSHED REF's copy of the workflow file. It is the only live CI host today.

Ruling, adopted: put the cross-host trio's cheap half on exact-gates now,
and give equivalence.yaml the same tokens for when it goes live.

- `exact-gates.yaml`: one step after the gates loop running
  `bcf-equivalence.R compare <baseline> --cross-host` and
  `multinomial-equivalence.R compare <baseline> --cross-host`. It sits
  OUTSIDE the `GATE_ARG` loop: that loop passes `quick` to every gate,
  and the baselines are recorded at FULL settings, which the settings
  guard ([[R/bartcore.R:458-472@736bfb05]]) refuses to compare against a quick run. Cost:
  measured on the box (section 4) and expected well under a minute per
  harness against the job's existing ~4-minute install and 30-minute
  timeout; the implementer substitutes the measured numbers in the step
  comment. Do NOT add `equivalence.R` here: its full 43-scenario x
  20-seed grid is what equivalence.yaml allots 120 minutes and
  `EQUIVALENCE_CORES` fork-parallelism to, and it would exhaust
  exact-gates' budget. The gaussian harness stays on equivalence.yaml.
- `equivalence.yaml`: add `--cross-host` at [[R/bartcore.R:88-89@736bfb05]] and [[R/bartcore.R:115-116@736bfb05]]. Zero
  runtime, zero new yaml, and the weekly cron starts firing at the
  merge - at which point equivalence.yaml carries the gaussian arm and
  exact-gates carries the per-push half.

So the standing cross-host gate after the box retires is exact-gates on
every push plus equivalence.yaml weekly from the merge. Until the slice's
own push, the box run is the whole gate.

## 6. Same-host regression safety

The default mode stays the 12-scenario identical-draws gate,
byte-for-byte. Four properties:

1. `--cross-host` is `setdiff`ed off `args` before `mode` is read
   ([[equivalence.R:34-35@736bfb05]] pattern). Absent, `crossHost` is FALSE and
   `exempt` is `character(0L)`.
2. `setdiff(names(a), c("summaries", character(0L)))` is value-identical
   to today's `setdiff(names(a), "summaries")` - same names, same order -
   so `channels`, `ok`, the `identical (all %d channels: %s)` line,
   `pointMismatch` and the terminal line are unchanged. Tier 1, tier 2
   and the `GATED` wording are all inside `if (crossHost)`.
3. The non-finite guard (2.4) runs unconditionally only under the flag;
   in the default mode it runs inside the mismatch branch, which the
   standing 12/12 bitwise gate never enters.
4. Derive-when-absent (section 3) is reachable only after a bitwise
   mismatch, so it changes output only on a compare that is already
   failing.

Slice gate, on the recording host (arm64/macOS): capture `compare`
stdout against `bcf-equivalence-6e3b9fb8.rds` before and after the patch
and `diff` - must be empty, 12 lines of `identical (all 7 channels: mu,
tau, glue, sigma, train, varcount, varcount.tau)` (8 channels on the five
verdict-carrying scenarios). Same for multinomial (11/11) and gaussian
(43/43 `--strict-coverage`). No baseline re-records at this slice.

## 7. RC-tip re-record: pending checklist

Mechanical, same-host, at the RC tip; nothing here can be done early.

1. Recording host only (arm64/macOS, the host every `current` MANIFEST
   row names). `R CMD INSTALL --preclean .` at the RC tip.
2. NEUTRALITY FIRST: `compare` against `bcf-equivalence-6e3b9fb8.rds`,
   expect 12/12 `identical`. If any scenario moves, STOP - the re-record
   then changes a recorded draw and owes a P17 ORACLE in its MANIFEST row
   (MANIFEST:16-23), not an adjudication.
3. `record benchmarks/baselines/bcf-equivalence-<rc-sha>.rds`, `<rc-sha>`
   the short hash of the ENGINE tip whose behavior the file records (the
   4d9a3337 precedent: the .rds names the commit a checkout can reproduce
   it from, not the records commit that follows).
4. SELF-REPRODUCTION: fresh `--preclean` install of that tree, `compare`
   against the new file, 12/12 identical.
5. CROSS-HOST CONFIRM: `--cross-host` against the new file from an x86
   leg - the box if it survives, otherwise the exact-gates run on the
   re-record push. Expect 12/12 `tier 1 PASS`.
6. FOUR-PLACE OBLIGATION, one commit (the 9a20bb4a precedent, which
   touched exactly the .rds, equivalence.yaml, MANIFEST, TODO and
   feature-matrix.md):
   - `benchmarks/baselines/MANIFEST`: new `current` row carrying the
     step-2 neutrality partition and the machine/arch string; demote the
     6e3b9fb8 row to `historical`.
   - `.github/workflows/equivalence.yaml`[[equivalence.R:89@23b9cde7]] AND the new exact-gates step
     from section 5 - the re-record now has TWO workflow pins, not one.
   - `TODO` ledger line.
   - FORWARD-facing docs pins only: `docs/design/feature-matrix.md`[[equivalence.R:647@23b9cde7]]
     (the `masked` scenario pin) and [f39] at [[equivalence.R:786-789@23b9cde7]]. Leave the
     historical landing records that name 6e3b9fb8 as a fact about a past
     tip - `docs/design/empty-leaf-veto.md`[[equivalence.R:329@23b9cde7]],
     `docs/design/multinomial-mutation-arc.md`[[equivalence.R:648@23b9cde7]] and [[equivalence.R:1066@23b9cde7]] - check each
     cite's tense before touching it.
7. ANCHOR SHIFT, easy to miss: inserting a MANIFEST row moves every line
   below it, and `docs/design` cites `MANIFEST:<n>` by line, checked by
   `tools/check-doc-freshness.R` (lint.yaml). The 6e3b9fb8 precedent
   repinned [f39] to `MANIFEST:16,42,49`; re-derive those numbers from
   the `git diff -U0` line map, not by a content pass. [f39] also carries
   per-baseline SCENARIO COUNTS - recompute them.
8. ALSO STALE and outside the four places:
   `docs/plans/review-2026-08-24/gate-ledger.md`[[equivalence.R:111@23b9cde7]] ("Current rows")
   still names `equivalence-5a3bc276.rds`. Fix it in the same commit.
9. Keep the old .rds as `historical`; nothing is deleted.

## 8. CI/consumer impact, residues, forks

CI impact, verified at this tip: `benchmarks/**` is `paths-ignore`d by
check-standard.yaml, cpp-tests.yaml, sanitizers.yaml and pkgdown.yaml;
exact-gates.yaml ignores `docs/**`, `TODO`, `**.md` and does NOT ignore
benchmarks, so a benchmarks-only push fires it. One correction to
"path-ignored except by exact-gates": lint.yaml ignores only
`benchmarks/baselines/**`, so `benchmarks/R/*.R` IS linted - the new code
must pass `air format --check .` and `lintr::lint()` against the slice's
own lib. Touching `.github/` fires all six push workflows, and a
benchmarks-touching push right after a slice push cancels the slice's
in-flight exact-gates run (`cancel-in-progress`): space the pushes.

Consumer impact: NONE. `benchmarks/` is `.Rbuildignore`d ([[equivalence.R:18@5a3bc276]]) and ships
nothing; no `R/`, `src/`, `inst/`, `man/`, `NAMESPACE` or
`inst/include/dbarts/dbarts.h` file is touched; no ABI, no hash re-bake,
no lockstep rebuild.

SCOPE, ruled: the flag lands in `multinomial-equivalence.R` in the SAME
slice, under the same two-tier verdict - `train`, `test`, `runVarcount`
gated ([[equivalence.R:103@5a3bc276]]), `forestFits`, `varcount`, `accepted`, `installed` exempt
([[equivalence.R:119-126@5a3bc276]] plus the call-site appends), same compare-loop shape
([[equivalence.R:595-646@5a3bc276]]), and its baseline already carries `summaries`. `equivalence.R`
needs NOTHING: its scenarios record a single seeds x summaries matrix and
its verdicts ride into that matrix (`recordVerdict`, [[equivalence.R:1057-1058@5a3bc276]], [[equivalence.R:1665@5a3bc276]]),
which is why it already exits 0 cross-host. `train`/`test` are
continuous, `runVarcount` combinatorial - the tier-1 split of section 2.2
applies unchanged.

Residues:

- [[benchmarks/README.md:55-56@5a3bc276]] says "Baselines are RNG- and
  build-dependent; regenerate them rather than reusing across machines",
  which `--cross-host` now contradicts. Rewrite those three lines to
  state the cross-host mode (a baseline IS reusable off-host under the
  flag, with the snapshot channels exempt and the two-tier verdict), and
  extend the [[benchmarks/README.md:59-66@5a3bc276]] statistical paragraph with the tier-2 weakness in one
  sentence.
- The Welch z at [[bcf-equivalence.R:508-513@5a3bc276]] (and
  [[multinomial-equivalence.R:600-605@5a3bc276]]) does `sa$mean - sb$mean` with no
  length check, so a channel SHAPE change silently recycles into a
  garbage z rather than erroring. That has already happened once - the
  `varcount` forest-axis widening (a16d703c). Add
  `stopifnot(length(sa$mean) == length(sb$mean))` inside the `Map`; two
  lines, fold it in here.
- `TODO`: close `equivalence-harness-statistical-mode` ([[multinomial-equivalence.R:71-90@a16d703c]]) at this
  landing - both REMAINING clauses are discharged - and reduce
  `bcf-baseline-cross-host` ([[multinomial-equivalence.R:33-36@a16d703c]]) to the RC-tip re-record. `rc-gate`
  ([[multinomial-equivalence.R:188-189@a16d703c]]) follows.
- The plan doc this becomes needs a `docs/plans/INDEX.md` row, or
  `tools/check-doc-freshness.R` fails lint.yaml.

Open fork (one, and it is genuinely open):

F1. `rtol` for tier 1. (a) Ship `1e-8` as designed and widen only if the
box run's observed ratio exceeds 0.01 - risks a second pass if the
observed deviation is larger than expected. (b) Set `rtol` from the box
measurement in the first place, leaving the constant blank until section
4 reports - costs a serialization point between the box run and the
final patch. RECOMMEND (a) with the calibration clause: 1e-8 is far
looser than a locked stream needs and far tighter than any defect class
survives, so the expected outcome is that no widening is needed and the
box number becomes the doc's evidence rather than its input.

## Landing note (2026-08-26)

LANDED at 3f532af2 (design record cbb1cb97). Both bcf-equivalence.R and
multinomial-equivalence.R gain a `--cross-host` compare-only flag under
the two-tier verdict: tier 1 LOCKED is the gate - a relative-deviation
bound (`rtol = 1e-8`, `atol = rtol * max|a|`) on continuous channels,
`identical()` on combinatorial channels; tier 2 `decoupled: statistical`
is a labelled-weak fallback, a Welch z with an ESS-adjusted denominator.
A non-finite guard runs on every gated channel; the partition assert
(section 1) fires both at script load and inside the compare loop over
`names(a)`. `exact-gates.yaml` runs both compares on every push,
outside its quick `GATE_ARG` loop; `equivalence.yaml` carried the same
`--cross-host` tokens for when it went live. Fixed-at 2518a0b5:
`equivalence.yaml`'s bcf/multinomial cross-host jobs were deleted as
duplicates of what exact-gates.yaml already runs on every push against
the same pinned baselines; only the gaussian statistical job and its
cron remain. `benchmarks/README.md` is rewritten at the "regenerate
rather than reusing across machines" sentence.

Gates: default mode byte-identical before/after the patch (bcf 12
identical, multinomial 11, equivalence.R 43 compared / 0 skipped);
same-host `--cross-host` 12 and 11 tier 1 PASS; lint parity (0 new
lints, 8 pre-existing brace_linter); both workflow yaml parse, `bash
-n` clean, no apostrophe in the new run block.

BOX validation (x86-64 Ryzen 7 3700X, R 4.3.3, gcc 13.3 - its final
engine leg before retirement): bcf 12/12 tier 1 PASS, worst ratio
3.7e-06; multinomial 11/11, worst 2.1e-06; equivalence.R 43 compared /
0 skipped, all max |z| 0.00; tinytest all ok 7478. `rtol` stayed at
1e-8 - the widening clause (section 2.2) did not fire, the observed
ratios sitting ~2700x under its trigger.

Discrimination probes, all run on the box. D-1' (`[[chain.hpp:5291@2518a0b5]]`, a
zero-RNG forest-0 `totalFits` write on a gated channel's own write
path) 12/12 tier 1 FAIL and 12/12 tier 2 FAIL, exit 1. D-2 (the same
mutant, unflagged) exit 1, naming `mu`/`tau`/`glue` in all 12 -
confirms the flag is load-bearing rather than a no-op. D-3
(`varcount` added to `snapshotChannels`) fires both asserts. D-4
(`[[R_interface_bartcore.cpp:2174@2518a0b5]]`, `sdModerate` x1.2 scoped to forest 1)
12/12 tier 1 FAIL, but tier 2 flagged only 11 of 12
(`set_predictor` max |z| 3.47) - the calibration evidence that tier 2
alone would pass a 20 percent node-prior widening. The independent
gate-runner reproduced the class on this host: a 1e-3 relative
perturbation of one channel fails tier 1 but tier 2 cannot distinguish
it, exit 0.

Recorded deviations: under the flag the three per-scenario lines
(exempt, non-finite, tier verdict) always print - the tier 1 line
appends "all N GATED channels bitwise identical" when so, rather than
a separate replacement pass line - so the path is exercisable
same-host without a real cross-host divergence. About 250 lines were
added per script against the ~90-line budget estimate, because the two
scripts are independent and duplicate the tier/guard helpers rather
than sharing them. The ubuntu-latest CI compare (`exact-gates.yaml`)
is now the standing cross-host gate, the box having retired.

Pending: the RC-tip same-host re-record of bcf-equivalence remains - a
refresh, not a prerequisite; checklist in section 7 of this doc,
including the MANIFEST anchor shift and `[[gate-ledger.md:111@2518a0b5]]`. Door: a
seeds-axis re-record, post-RC.
