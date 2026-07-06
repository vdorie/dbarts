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
