# Installing a stored state: the checked form and the forced form

Status: LANDED 2026-10-10: Part A (missingness first seen, 2269918e) and Part B (setState's two forms, 062ef59a).
Plan: [setstate-force-update.md](../plans/setstate-force-update.md). Rulings: dec-B305, dec-B310, dec-B318,
dec-B398 to dec-B402 and dec-B418 in docs/decisions.md.

A state is the chain, and installing one puts a chain into a sampler that has its own model and its own data.
Most states fit. Some hold a tree the sampler cannot hold as it is stored. This note says what `setState` does
with each, and is the home of those rules; what a state carries is
[state-not-model.md](state-not-model.md)'s, and which columns can hold a missing value is
[mia-missingness.md](mia-missingness.md)'s.

## The two forms

`setState(newState, forceUpdate = FALSE)` asks first. A state every tree of which fits is installed and the
call returns TRUE. A state with a tree that would have to be changed is not installed: the call returns FALSE
and the sampler is as it was, its trees, drawn values, latents, generators, cut grid, kept draws and cached
`state` field, and a pointer that was dead stays dead. Both values are returned visibly, as an unforced
`setPredictor`'s are.

`setState(newState, forceUpdate = TRUE)` installs, repairing what does not fit, and returns NULL invisibly
whether or not anything was repaired (dec-B310: the unforced value says whether a statistically valid update
was made, which has no meaning once the update is forced). `copy()` and a reload force and report nothing,
having no other state to fall back on; so does stan4bart's reload (dec-B402).

A state that is broken, or of a shape that does not fit, is an error in both forms: not a state, an older
encoding, another count of chains, forests or trees, another leaf model, a malformed block or tree, a cut grid
that repeats a point, latents the family cannot hold, saved gp draws under other lengthscales, saved-tree
blocks that name different numbers of draws.

## What does not fit, and its repair

Only the live trees are judged, mean and variance, each built from its stored form and partitioned over the
sampler's rows. A kept draw is replayed over new rows and never held to this partition.

| a tree that | forced |
|---|---|
| has a leaf no row reaches | the split above it becomes one leaf with everything beneath |
| holds a split outside the range the splits above it leave | the same |
| records a side for missing values on a column, not pooled, that has never held one here | the side is dropped |
| splits on a column its forest may not use (`vars`, a `blocks` group, a moderator's columns, a restricted variance forest) | cut back at the first such split from the root |
| breaks its forest's `interactions` limit | cut back at the first such split from the root |
| is monotone with leaf values out of order as stored | every leaf set to 0 |

A merged leaf takes the mean of the leaves it replaces, weighted by their rows under the state's own working
weights, geometric in the variance forest. On each path from the root the first offending split goes and the
splits above it stay. A merge can leave a monotone tree out of order, and it is then reseeded as well. The
last three rows are repaired under dec-B398 ("Repair, as we do elsewhere."). A warm start is not a state
install: it refuses a donor that breaks a column restriction or an interaction limit, and reseeds one out of
the monotone order.

The verdict and the install stage a tree by one routine,
[`Tree::stageFromFlat`](../../src/bartcore/tree.hpp), the verdict on a scratch tree
([`Chain::stateIsValid`](../../src/bartcore/chain.hpp)) and the install on the live one
([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp),
[`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp)), so the unforced call declines exactly the
states a forced install reports as altered. The repair is one pass of
[`Tree::collapseEmptyNodes`](../../src/bartcore/tree.hpp), whose walk also takes a forbidden split
([`Tree::splitIsForbidden`](../../src/bartcore/tree.hpp)).

The verdict builds and partitions every tree once more than a forced install does, so an unforced restore
costs more. Measured with 200 trees on 10 columns, one chain and one thread (bench-sampler.R's `restore`
grid, 2026-10-10, on a machine under load, two runs within one percent), an unforced restore took 1.64 times
a forced one at n = 1000 and 1.72 times at n = 1e5, which is 0.9 and 1.3 sweeps more; a forced restore cost
what a restore cost before the check existed (1.000 and 1.001 times it). A loop that restores its own state
onto rows it has put back, where nothing can need repair, passes `forceUpdate = TRUE` and pays that.

## What is clean whatever it differs in

- The response mapping. Nothing is converted on install (dec-B418), so a state stored before a re-anchor or by
  a sampler on another response goes in as stored and is clean.
- The latents. Under other case weights or another censoring status they are redrawn after the install, from
  the installed generators.
- `k`. It is installed as recorded whatever leaf prior the recipient holds (dec-B401); the leaf sd is then
  `k.scale / k` against the recipient's `k.scale`. A recipient holding `k` fixed keeps its own.
- The kept draws (dec-B318: they are output, not model or chain). A state from a store the size of the
  sampler's has its ring copied slot for slot, so a copy continues its source's layout. From a store of any
  other size the newest recorded draws that fit go to the sampler's first slots, oldest first; `n.samples` and
  the store's size do not change, and a sampler keeping none takes none
  ([`Sampler::setState`](../../src/bartcore/sampler.hpp)). A re-created sampler still takes the state's
  capacity, having no other record of it.
- A value the sampler holds fixed, DART weights, and a generator of another kind, each left as the sampler's.

## Not in this note

A state stored before a column's first missing value and installed after it draws its directions (dec-B435);
that is its own slice. `installTrees`'s `forceUpdate`, and its NULL, come with the warm start's salvage.
