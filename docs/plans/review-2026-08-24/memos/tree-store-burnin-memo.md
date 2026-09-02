# keepTrees store: predict is rotated across successive run() calls

Scope memo, tip 04838898; probes from a private copy of it. `git diff fab82735
04838898` touches only TODO, a docs/plans file and tests/cpp/test_fuzz.cpp, so
the fab82735 report carries over. Scripts:
`<scratchpad>/burnin-scope-32065/{repro,repro2..5,capi-probe}.R`.

## 0. Correction to the report

Headline right (predict's draws rotate, content correct, silently); mechanism
wrong twice.

- **Burn-in sweeps are NOT recorded.** `[[chain.hpp:1365-1367@fab82735]]` gates on
  `record = (iteration+1) % numThin == 0 && iteration/numThin >= numBurnIn`, and
  every store site is `if (record && ...)` (`[[chain.hpp:1494-1498@fab82735]]`, `[[chain.hpp:1537@fab82735]]`,
  `[[combiner.hpp:1375-1382@fab82735]]`); `[[sampler.hpp:427-429@fab82735]]` advances the cursor by
  `numSamples` only; `[[tests/cpp/test_state.cpp:210@fab82735]]` already pins it.
- **`sampler$run(n.burn, n.samples)` on a fresh sampler is CORRECT**, as is
  `dbartsControl(keepTrees=TRUE, n.burn=2, n.samples=6)` + `$run()`.

Real trigger: the sampler write cursor `currentSampleNum_` carries across
`run()` calls while all four readers index slots 0..capacity-1. Rotation = the
cursor at the sampling run's start = (draws recorded before it) mod capacity.
"Rotate by n.burn" is the case where the host ran burn-in as its own RECORDED
run `run(0L, n.burn)` with keepTrees left on - the natural hand-driven idiom and
the shape of `runWithBurnIn` (`[[R/bart.R:550-590@fab82735]]`) minus its keepTrees toggle.
Fires at `n.burn == 0` too.

## 1. Reproduction

Seeded single-forest gaussian, n=60 p=6 n.trees=5 n.chains=1 n.threads=1,
`keepTrees=TRUE`; statistic `max|s$predict(x.train) - run()$train|`, then k.

| probe | form | k=0 | best k | resid |
|---|---|---|---|---|
| a | fresh `run(2L,6L)` / all-in-control `$run()` / 2nd+3rd `run(2L,6L)` / `run(7L,4L)` cap 4 | 7.1e-15 | 0 | 7.1e-15 |
| **A** | `run(0,nb)` then `run(0,6)`, cap 6, nb=1..4 | 6.65/6.57/6.07/6.12 | **1/2/3/4** | 7.1e-15 |
| **B** | `run(0,7)` then `run(0,4)`, cap 4 (nb > ns) | 4.719 | **3** = 7 mod 4 | 7.1e-15 |
| **C** | `run(2,4)` x3, cap 6 (ns != cap) | - | **0, 4, 2** | 7.1e-15 |
| **D** | 2 chains, `run(0,2)` then `run(0,6)` | 6.57 / 10.70 | **2 and 2** | 7.1e-15 |
| E | A(nb=2) then `storeState`/`setState`, or `$copy()` | 6.567 | 2 | `max|P1-P0| = 0` |
| G | A(nb=2), keepTrees toggled off/on between | 7.1e-15 | 0 | - |
| **capi** | flat C API `capi_run(0,3)` then `capi_run(0,7)`, cap 7 (single `capi_run(3,7)`: 0) | **0.5204** | **3** | **0** exact |

Slot map pinned directly (A, nb=2): `max|P[,1]-train[,5]| = 7.1e-15`,
`max|P[,3]-train[,1]| = 3.6e-15` - output slot i holds draw `(i-base) mod cap`.
Rotation is per-SAMPLER, not per-chain (D): one base fans to every chain
(`[[chain.hpp:2582-2584@fab82735]]`). It survives serialization (E): `state.currentSampleNum`
round-trips (`[[sampler.hpp:666@fab82735]]`, `[[sampler.hpp:731@fab82735]]`), readers ignore it. Partial fill (cap 6,
one `run(0,3)`) gives 6 draws, 3 recorded plus 3 zero-leaf - a separate edge.

**Live consumer that DEPENDS on the carry**: `bart2(family="nbinom")` loops
`bartcoreRun(bc, if (s == 1L) control@n.burn else 0L, 1L)` per kept draw
(`[[R/bart.R:2068-2070@fab82735]]`); correct today, `predict(fit,x,"bart")` vs
`extract(fit,"bart")` best k=0 resid 1.9e-15. Ordinal (`[[R/bart.R:1820@fab82735]]`) runs once.

## 2. Mechanism

**Store and write.** Per-forest circular buffer, slot-major
`capacity x numTrees` (`[[chain.hpp:426-432@fab82735]]`, `[[combiner.hpp:221-232@fab82735]]`), capacity
fixed at creation from `control@n.samples` (`[[R_interface_bartcore.cpp:1866@fab82735]]`) or
by `setTreeStorage` (`[[R_interface_bartcore.cpp:4881@fab82735]]`), allocated whole up front (`[[chain.hpp:2546-2580@fab82735]]`).
`Sampler::run` sets each chain's base to the sampler cursor
(`[[sampler.hpp:272@fab82735]]`); the chain writes slot `(savedSlotBase + sampleNum) %
capacity` (`[[chain.hpp:1494-1498@fab82735]]`; variance twin `[[chain.hpp:4244-4248@fab82735]]`; the amplitude
rescale patches the same slot in place, `[[combiner.hpp:1375-1382@fab82735]]`, so the glue is
folded INTO the saved leaves and cannot itself mis-pair), and the cursor then
advances by `numSamples` (`[[sampler.hpp:427-429@fab82735]]`). Reset only by `setTreeStorage`
when keepTrees/capacity change (`[[sampler.hpp:958-970@fab82735]]`) and by a warm-start install (`[[sampler.hpp:913@fab82735]]`).

**Read (the defect).** All four readers walk slots 0..capacity-1 and never
consult the base: `[[sampler.hpp:548@fab82735]]` (predictColumns), `[[sampler.hpp:601@fab82735]]`
(predictPerForestColumns), `[[sampler.hpp:650@fab82735]]` (predictVarianceColumns),
`[[R_interface_bartcore.cpp:7275-7281@fab82735]]` (getTrees; its 0-based `sampleIndices` are
range-checked against capacity at `[[R_interface_bartcore.cpp:7475@fab82735]]`, `[[C_interface.cpp:842@fab82735]]`). The bridge
sizes predict's draw axis as `capacity`, not the last run's `numSamples`
(`[[R_interface_bartcore.cpp:5615-5617@fab82735]]`, `[[R_interface_bartcore.cpp:5783@fab82735]]`).

**bart()/bart2()/rbart_vi dodge.** `[[R/bart.R:1166-1171@fab82735]]` forces keepTrees FALSE
whenever `n.burn > 0`; `runWithBurnIn` (`[[R/bart.R:550-590@fab82735]]`) runs burn as
`sampler$run(0L, n.burn)` with storage OFF (capacity 0, so `[[sampler.hpp:428@fab82735]]`
skips the advance), then re-installs keepTrees via `setControl` - flipping
capacity 0 -> n.samples resets the cursor (`[[sampler.hpp:969@fab82735]]`) - and runs from
base 0. `rbart_vi` toggles identically (`[[R/rbart.R:939-941@fab82735]]`, `[[R/rbart.R:969-972@fab82735]]`). Probe
G confirms the toggle is what saves them.

**Docs.** `[[docs/design/core-generalization.md:680-683@fab82735]]` specifies the store
("capacity = n.samples at creation ... currentSampleNum advancing per run") but
is silent on read order; `[[docs/design/multinomial-mutation-arc.md:586-589@fab82735]]`
records the carry as deliberate and measured;
`[[docs/plans/correctness-audit.md:433@fab82735]]` asserts keepTrees indexing "is uniform
across store, read, and predict including circular wrap" - the false claim.

## 3. Options

**(A) Don't record burn-in.** Already true, pinned by `[[test_state.cpp:210@fab82735]]`;
zero lines, fixes nothing. **(D) Refuse keepTrees with n.burn > 0 on the
hand-driven path** does not address the bug either - probe (a) shows `run(2,6)`
already correct, probe (C) shows the rotation at `n.burn == 0` - and would
refuse a legitimate call. Both dropped.

**(B) Readers iterate from the cursor. RECOMMENDED.** Output index i reads slot
`(currentSampleNum + i) % capacity`, so slot order == chronological order over
the `capacity` most recent recorded draws. Files: `[[sampler.hpp:548@fab82735]]`, `[[sampler.hpp:601@fab82735]]`,
`[[sampler.hpp:650@fab82735]]` (one expression each plus a `savedSlotForDraw(i)` helper beside
`currentSampleNum()` at `[[sampler.hpp:476@fab82735]]`), and `[[R_interface_bartcore.cpp:7275-7281@fab82735]]`.
`SamplerBase::currentSampleNum()` already exists as a virtual (`[[facade.hpp:202@fab82735]]`,
`[[facade.hpp:480@fab82735]]`): **no new virtual, no facade churn, no dbarts.h signature change**, no
`stateFormatVersion` bump (stays 2, `[[facade.hpp:6250@fab82735]]`). Preserves the nbinom loop and the
Door 3 batching equivalence (both end with cursor 0, so their read order is
unchanged). Cost: `[[tests/cpp/test_state.cpp:225-240@fab82735]]` asserts a partial run's new
draws land at output 0..1; under B they land at the tail (~6 lines).
`[[dbarts.h:837-849@fab82735]]` and `?dbartsSampler` each gain one prose sentence.

**(C) Reset the cursor at each sampling phase.** One line at `[[sampler.hpp:272@fab82735]]`.
**Disqualified**: breaks `bart2(family="nbinom")` outright (`[[R/bart.R:2068-2070@fab82735]]`
records draw s in its own `run` call; with a reset every draw lands in slot 0),
deletes the accumulation idiom the ring exists for, contradicts `[[R/bart.R:231@fab82735]]` above.

**(B') B plus a fill counter**, so predict reports only the draws actually
recorded (also fixing partial fill). Adds a `size_t` to `SamplerState` (format
2 -> 3) and changes predict's returned extent - a visible surface change; defer,
recording it as the open edge in the design doc.

**Equivalence trio: BITWISE under B, no re-record.** `bcf-equivalence.R` and
`multinomial-equivalence.R` contain no `keepTrees`/`predict`/`$run` at all. In
`equivalence.R` the hand-driven scenarios use `sampler$run(nskip, 0L)`
(numSamples 0 - no advance, `[[sampler.hpp:428@fab82735]]`) then one recorded run
(`[[sampler.hpp:1034-1087@fab82735]]`), summarizing only run channels, never predict; the two that DO
read the store, `fitViaHazard` (`[[sampler.hpp:1186@fab82735]]`) and `fitViaHurdle` (`[[sampler.hpp:1224@fab82735]]`), go
through `bart2` with `n.burn > 0`, i.e. cursor 0 and `numSamples == capacity`,
so `(0+i) % cap == i`. Holds for equivalence-5a3bc276.rds (43 strict-coverage),
bcf-equivalence-6e3b9fb8.rds (12), multinomial-equivalence-4d9a3337.rds (11).

## 4. Tests

New `inst/tinytest/test-sampler-keeptrees-order.R` plus additions to
`tests/cpp/test_state.cpp`. Identity throughout, to 1e-12:

1. R5 two-phase `run(0L,nb)` then `run(0L,ns)`, nb in 1:4 with ns == capacity
   and nb > ns; plus `dbarts()` + `$run(nb, ns)` in one call (guard on the
   currently-correct path); plus `n.chains = 2`, per chain slab.
2. Two successive recorded runs, `ns != capacity` (`run(2,4)` x3, cap 6): pin
   the whole store chronologically, not just the last run.
3. `storeState()`/`setState()` round trip of a store written by (1): predict
   identical before and after, still equal to `$train`.
4. `getTrees(sampleNums = s)` agrees with the predict draw at index s; nbinom
   `predict` vs `extract` - the guard against option C creeping back.
5. C++: after `run(0,2); run(0,n)` on a capacity-n store, `predict` equals the
   second run's recorded `testFits` exactly; update `[[test_state.cpp:225-240@658869ac]]`.
6. Planted mutation: revert `savedSlotForDraw` to identity; (1), (2), (5) must
   fail at O(response sd), only the rotated form passing.

## 5. Price

Raw added lines, 1.5x stop each: `sampler.hpp` 12; `R_interface_bartcore.cpp`
8; `dbarts.h` 3 (prose); `docs/design/core-generalization.md` 10;
`tests/cpp/test_state.cpp` 25; `inst/tinytest/test-sampler-keeptrees-order.R`
90; `man/dbartsSampler.Rd` 4. **Total ~152, stop 228.** Gates: `tests/cpp` full from clean plus `./test_sampler`; full tinytest FAILURES == 0
after `R CMD INSTALL --preclean` into a private lib; the trio with the bitwise
verdict above; ASAN/UBSAN (engine code moves in a path predict reaches);
`air format --check` + `lintr`; `R CMD check --as-cran` from a clean-copy
tarball. `bench-sampler` NOT owed - predict is off the sampling path.

## 6. Risks

- **Semantics call, not just a fix.** B redefines predict's draw axis as "the
  `capacity` most recent recorded draws, oldest first" - public contract, must
  land in `dbarts.h` and `?dbartsSampler`, unambiguous only once partial fill is
  answered (B'). "Last run's draws only" is a different design, format bump and
  changed return extent. `[[docs/plans/correctness-audit.md:433@658869ac]]` is false and
  should be corrected in the same change.
- **stan4bart** is on the bartcore branch and calls `dbarts_sampler_predict`
  after its own run loop; if it drives more than one recorded run per fit it
  reads a rotated store today - one grep owed there before landing.
- Partial fill stays sharp under B: predict still returns `capacity` draws
  including zero-leaf slots, at the head rather than the tail.
- Measured on gaussian constant-leaf single-forest plus the nbinom loop;
  BCF/multinomial/heteroscedastic share the same three slot loops, so the
  analysis carries but was not measured.
