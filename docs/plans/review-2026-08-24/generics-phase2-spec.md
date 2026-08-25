# Generics slice, phase 2 spec (orchestrator ruling on the survey's ambiguities)

Survey: j2-survey.md sections 0 and 4. Ambiguities in section 5 DECIDED as follows.

SHARED RULES. plot: left panel(s) = trace of the family's scalar posterior parameter(s),
right = observed-vs-fitted posterior-interval panel with a reference line; keep plot.bart's
plquants / cols formals and by-name refusals; par(mfrow) set and restored like plot.bart.
loglik: extract-only; sample = "test" refused by name; NO weight term (all four families
refuse weights at fit time); NO offset term (offsets sit inside the stored channels);
combined = S x N, uncombined = chains-first; zero-weight NaN rule as the gaussian arm.
A combineChains formal (default TRUE) on extract.bartMultinomial / .bartOrdinal / .bartNegbin,
every channel routed through reshapeChainedChannel / combineOrUncombineChains.
as_draws: IMPLEMENTED for all four (overrides plan D6): expose the fit's scalar posterior
parameters (what summary() tables), never per-observation channels; registered the way
R/hooks.R registers bart's.

bartNegbin. plot mfrow(1,2): P1 dispersion trace (plotSigmaTrace body, ylab "dispersion (r)",
no burn-in segment - documented - type = "s"); P2 observed y vs posterior median of
yhat.train with plquants bars, abline(0,1), ylab "posterior interval for E(Y|x)".
loglik: dnbinom(y_i, size = dispersion[s], mu = yhat.train[s,i], log = TRUE), dispersion
paired chain-fastest like sigma. dims = dim(ev). as_draws vars = c("dispersion","sigma","k","tau").

bartOrdinal. plot mfrow(1,2): P1 one trace per FREE cutpoint (gamma_2..gamma_{K-1}); at
K = 2 skip P1, P2 becomes the single panel. P2 latent.train median + plquants bars, obs
ordered by median eta, coloured by observed level, dashed hlines at posterior-median
cutpoints. loglik: log(yhat.train[s,i,k]), k = as.integer(y_i) (stored form). dims = ev
without the K margin (= ppd shape). as_draws vars = c("cutpoints","sigma","k","tau"),
cutpoint[1] kept (pinned at 0, documented).

bartMultinomial. plot mfrow(1,2): P1 existing per-category trace; P2 label ingestion =
median and plquants interval of p(observed category) vs the median, abline(0,1), xlab
"median of p(observed category)"; counts ingestion with multi-trial rows
(is.matrix(y) && any(rowSums(y) > 1)) = gaussian-panel form over n*K cells (observed
proportion y_ik/n_i vs interval of p[s,i,k], identity). loglik: lgamma(n_i+1) -
sum_k lgamma(y_ik+1) + sum_k y_ik log p[s,i,k] (= dmultinom, coefficient INCLUDED), one
formula for both ingestions; unit = the observation row (Rd: loo is leave-one-row-out);
dims = ev without K. Comment that the engine's NaN per-observation flag for this family stays
by design. as_draws vars = "meanProb" (K variables meanProb[<level>]).

bartHurdle. plot mfrow(2,2): P1 plotSigmaTrace(x$positive$first.sigma, x$positive$sigma)
(helper's own mfrow made optional); P2 occupancy median pi vs interval, ylab "posterior
interval for P(Y > 0 | x)"; P3 observed log(y) on y > 0 rows vs interval of f there,
identity; P4 observed y over ALL n vs interval of E[y|x] = pi * exp(f + sigma^2/2), identity.
loglik via hurdleParts(): y = 0: log1p(-pi); y > 0: log(pi) + dnorm(log y, f, sigma,
log = TRUE) - log(y). NO truncation. Natural-y scale WITH the Jacobian (DECIDED). Rd states:
(a) differs from a log-y computation by -sum_{y>0} log y_i, comparable to other models OF y;
(b) NOT the sum of the components' own loglik channels; (c) pi = P(Y > 0), complement of
brms's hu; (d) a count hurdle (door) would need truncation, this one must not. dims = S x n.
as_draws: union of components' scalar fields with "occupancy." / "positive." prefixes
(dot, not bracket), vars = c("sigma","k","tau").

Rd: man/bart.Rd's loglik dim contract restated as "the ev shape without the trailing category
margin" where it applies; test-pointwise-loglik.R dims adjusted. Family items in bart2.Rd's
own-class paragraphs (or wherever plot.bartMultinomial is documented).

PINS (red-then-green; mutation proofs): per family loglik == hand formula to 1e-12
(dnbinom / dmultinom incl. coefficient / log Phi-difference from latent + cutpoints /
log1p(-pi) + dlnorm composition); combined vs uncombined shapes; plot runs on pdf(NULL) for
every family incl. K = 2 ordinal and multi-trial-count multinomial; as_draws variable names;
refusals by name. Mutations: drop multinomial coefficient; invert hurdle pi; drop the Jacobian;
swap size/mu in dnbinom.

Doors (landing note only, NOT code): negbin rootogram panel; ordinal log_diff_exp tail
precision; negbin burn-in dispersion channel.
