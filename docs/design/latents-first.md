# Latents first: a sweep draws its augmentation variables before the trees

Status: PLANNED 2026-11-05, on the latents-first branch, not landed ([latents-first: a sweep draws its augmentation variables before the trees](../plans/latents-first.md#latents-first-a-sweep-draws-its-augmentation-variables-before-the-trees)); the user-facing text is the embedding page's "When the augmentation variables are drawn".

The class: POSTERIOR-CHANGING for an embedded sampler whose outer step integrates the latents out, and for a logistic sampler handed a new offset; draw-changing and posterior-preserving for a plain fit of a family with augmentation variables; bit for bit for gaussian, heteroscedastic gaussian, multinomial and an aft fit with no censored row.

## The question

Where in a sweep are the augmentation variables drawn - the probit and ordinal latents, the Polya-Gamma precisions of logistic, binomial and negative binomial, the Student-t scale mixing, the aft imputed times - so that a change the code driving the sampler makes between sweeps reaches the next sweep's trees?

Drawn after the trees, the next sweep's trees fit values drawn under the state before the change. The change is an offset, a predictor, sigma or a basis written between two `run(0, 1)` calls. What that does to the posterior depends on the outer step.

## The two outer steps

An outer step that CONDITIONS on the latent runs the sampler, reads the latent with `getLatents`, and draws its own quantity given that value (stan4bart's, `rbart_vi`). It is exact in either order: the latent it reads is the one the trees were drawn from.

An outer step that INTEGRATES the latent out scores a proposal with the marginal likelihood (the probit likelihood `pnorm(f)` of the response, say), installs what it accepts, and runs the sampler. The next sweep must fit trees to a latent drawn under the accepted value. Under the old order it fits one drawn before the change, so the sampler's conditional is not the one the outer step assumed, and the joint chain is off: on a small probit model with the latent averaged out, about 15 standard errors in a gate of 400,000 sweeps. Drawn first, the latent of the next sweep is a draw from its conditional given the fit and the inputs the trees will see, and the chain is exact.

Drawing first is exact under either outer step, whatever a setter changed, with no extra draw a sweep. The setters that drew to patch the old order (`setResponse`, `setActiveRows`, logistic `setWeights`) no longer do: each stores what it is given and the next sweep draws.

## The order

Per sweep: the per-sweep hook (where a caller makes its changes); the augmentation draw, against the combined fits as the previous sweep left them, each family's pairs in the order it always had (cutpoints then z, shape then omega, lambda then nu); the working response and weights pointers read again; mean weights under a variance forest; the level shift; the probit rescaling step; the trees; sigma; the combiner's glue; the variance forest; k and DART; the record and the draw callback. The step list is [One sweep](../architecture.md#one-sweep).

The draw precedes the rescaling step, not follows it. The rescaling step scales z, the leaves and k by one factor drawn from its conditional given them, so it reads z. A z drawn after it was drawn before the caller's change, and the rescaling step would then act on a stale latent: a build with the draw placed after the rescaling step is as far off as drawing last, which the gate measures.

Reading the working pointers every sweep, in place of once a run, is what makes a setter in the hook take effect in the same sweep.

## What a getLatents reader sees

After a run, the values the last sweep's trees were drawn from; they and the fits read back are a joint draw. Between a setter and the next run, the values from before the setter, except that a new data object starts them as creation does (an ordinal or negative binomial sampler from the kept cutpoints or shape) and an aft response leaves a censored row at its new bound. In the per-sweep hook, the previous sweep's. A callback that reads them on an inline run sees this sweep's; `dbarts_draw` carries none.

Under this order no next sweep reads the previous latent (each family's conditional draws from the fits, sigma and its family parameters, and the collapsed draws integrate the latent out), so a stored latent serves `getLatents` after a restore and nothing else, and a restored or copied sampler continues the original bit for bit. The state therefore no longer records which rows were at weight zero, nor fingerprints of the weights and the censoring: their only use was redrawing at a restore. A weight change is undone by `setWeights` and `setState` in either order, the sampler then continuing bit for bit as one that installed the stored state under the old weights would (not the uninterrupted chain, which a restore matches to about 1e-14).

A sampler keeps the ordinal cutpoints, the Student-t degrees of freedom and the negative binomial shape through `setData`, as `setResponse` keeps the cutpoints: the next sweep draws the latents from the kept values, a valid state under any data, and a caller who swaps data every sweep no longer holds those parameters at their starting values.

## The gate

[`outer-step-exact.R`](../../benchmarks/R/outer-step-exact.R) runs the real engine on 32 rows in four cells with one tree, enumerates the tree space and computes the exact posterior of an offset coefficient and of each cell's fit. An outer step draws the coefficient every sweep, in four arms: probit and logistic, each with the outer step integrating the latent out and with it conditioning on the latent. Pass: every cell fit and every grid point of appreciable mass within 4.5 standard errors. The gate checks where the draw sits, before or after the rescaling step, and nothing else; a stale working pointer or a missed invalidation is the unit tests' job.

## Consumers

stan4bart and `rbart_vi` condition on the latent: their draws move and they stay exact. A caller that integrates the latent out (treatSens's treatment model, the ability step of bairrtt, the vignette's fourth recipe) is exact where it was off. bartCause's probit and BCF fits draw differently at the same seed.
