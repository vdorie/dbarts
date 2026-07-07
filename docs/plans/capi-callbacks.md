# capi-callbacks

agent: opus
rng: neutral
budget: ~300 lines across dbarts.h, C_interface.cpp, facade, consumer test

## Goal

A host embedding BART in an outer Gibbs sampler can register a
conditioning callback invoked once per sweep inside a multi-iteration
run(), instead of round-tripping setSigma/setOffset + run(0, 1) per
outer iteration. An observer hook and the grouped results fields ride
the same additive API bump.

## Context

- The per-sweep dispatch point already exists: Sampler::run's
  pollInterrupt std::function (src/bartcore/sampler.hpp:180-186);
  dbarts_sampler_run does not forward it (src/C_interface.cpp:58-72).
- Worker threads must never call into R (sampler.hpp:253-278; the
  ProgressSink queue exists for this reason).
- dbarts_results lacks tau/groupEffects, which bartcore::Results
  carries (src/bartcore/chain.hpp); public-surface.md section 6
  deferred callbacks "until a consumer exists".
- Entry points are name-looked-up CCallables, so evolution is additive;
  bump DBARTS_C_API_VERSION.

## Decision (memo 2026-07-07; VD sign-off gates implementation)

1. dbarts_results extension. The struct is caller-owned, so growing it
   is safe only while nothing is released: stan4bart's
   IterableBartResults (bart_util.hpp:27) holds it as an uninitialized
   member and assigns the six fields by name - new fields would be
   garbage pointers without a one-line fix (value-initialize
   `current`). Post-release, growth is never safe: old-compiled
   callers zero-init only the old sizeof, so the new library would
   read stack garbage past the end. Recommendation: extend the struct
   with tau/groupEffects NOW (both packages unreleased, lockstep
   recompile; patch stan4bart's one line on its dbarts-1.0 branch),
   and freeze the struct at release - post-1.0 additions use new
   entry points or registration setters instead. Alternative: freeze
   at six fields today and add
   dbarts_sampler_setGroupedResultBuffers(sampler, tau, groupEffects)
   - no stan4bart edit, but the grouped outputs live off-struct
   forever.
2. Callback contract. Signature int cb(void* userData,
   dbarts_sampler*, size_t chainIndex, size_t sweepIndex, int
   isBurnIn), return 0 to stop; invoked on the calling thread BEFORE
   every sweep (unthrottled - pollInterrupt's 100ms throttle stays
   separate), which reproduces the host's setState-then-run(0, 1)
   loop exactly and serves the observer role (it sees the previous
   sweep's state; final state is read after run returns).
   Registration refused when chains run on worker threads. Caveat to
   sign off: inline multi-chain runs chains SEQUENTIALLY (chain 0
   completes all sweeps first) - per-chain RNG keeps draws identical
   to the interleaved loop, but a callback cannot condition chain c
   on chain c+1's progress; stan4bart's shape (one chain per sampler)
   is unaffected. Recommendation: allow inline multi-chain with the
   chainIndex argument and document the order; restricting v1 to
   numChains == 1 is the conservative alternative.
3. The rbart_vi custom-prior loop fix rides step 1's facade hook
   (bridge-internal, not the public C API); its follow-up plan is
   written when this lands.

Resolved (VD, 2026-07-07): extend the struct now, with the
release-time freeze (post-1.0 additions use new entry points);
callback contract as recommended - inline multi-chain allowed with
the chainIndex argument and the sequential order documented.
Consequences: DBARTS_C_API_VERSION stays 1 (nothing has shipped;
the extension is part of the initial contract, so step 2's "bump"
is void), and stan4bart takes the one-line value-initialization fix
(bart_util.hpp, `dbarts_results current;`) in lockstep before
submission.

## Constraints

- Callback contract v1: invoked on the calling thread, between sweeps,
  only when chains run inline (numChains == 1 or numThreads == 1) -
  stan4bart's shape. Multi-threaded chains refuse registration with a
  clear error; lifting that needs a barrier design, out of scope.
- Callback signature must let the host mutate conditioning state
  (sigma, offset) through the existing setters from inside the callback
  and signal early stop via return value.
- Out of scope: R-level exposure; replacing the run(0, 1) pattern in
  stan4bart itself.

## Steps

1. Facade: add a per-sweep hook alongside pollInterrupt (single
   std::function<bool(size_t iteration)> is enough for both roles;
   false = stop).
2. dbarts.h: dbarts_sampler_setCallback(sampler, fn, userData) with the
   threading contract documented in the header; extend dbarts_results
   with tau/groupEffects (NULL skips, matching the existing fields);
   bump DBARTS_C_API_VERSION.
3. C_interface.cpp: register the CCallable; wire the results fields.
4. Extend inst/tinytest/capi/consumer.c: drive a run with a callback
   that sets sigma per sweep; assert equality with the equivalent
   setSigma + run(0, 1) loop; assert early stop; read tau on a grouped
   fit if cheaply constructible, else skip.

## Verification

- test-capi.R passes (self-gates on toolchain).
- Full tinytest suite; component tests.
- Callback path vs setter-loop path: identical draws (the equality
  assertion in the consumer test is the gate).
