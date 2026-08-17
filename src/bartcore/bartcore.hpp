#ifndef BARTCORE_BARTCORE_HPP
#define BARTCORE_BARTCORE_HPP

// Header-only conjugate BART MCMC engine: constant, monotone, linear, and
// Gaussian-process leaves; ordinal and categorical split rules; Gaussian,
// probit, logistic, AFT (survival), ordinal, and negative-binomial response
// families, with Student-t continuous errors by scale mixture; heteroscedastic
// variance forests; grouped random intercepts; DART split-variable selection;
// BCF two-forest and multinomial softmax forest couplings; and sparse-column
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
