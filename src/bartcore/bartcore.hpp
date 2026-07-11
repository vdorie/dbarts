#ifndef BARTCORE_BARTCORE_HPP
#define BARTCORE_BARTCORE_HPP

// Header-only conjugate BART MCMC engine: constant, linear, and
// Gaussian-process leaves; ordinal and categorical split rules; Gaussian,
// probit, and logistic response families; grouped random intercepts; DART
// split-variable selection; BCF two-forest samplers; and sparse-column
// ingestion. See docs/design/core-generalization.md.

// IWYU pragma: begin_exports
#include "data.hpp"
#include "tree.hpp"
#include "model.hpp"
#include "scan.hpp"
#include "grow.hpp"
#include "moves.hpp"
#include "chain.hpp"
#include "sampler.hpp"
#include "facade.hpp"
// IWYU pragma: end_exports

#endif  // BARTCORE_BARTCORE_HPP
