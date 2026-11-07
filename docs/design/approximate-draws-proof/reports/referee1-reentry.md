I checked all four points myself: by reading, by an exact toy, and against the engine source on main. Points 1 to 3 hold, with three corrections. Point 4 is confirmed, but setData's limit is wider than report5 says. Toy: the referee's own scripts (not kept) Its state is theta in {1,2}, active in {0,1}, omega in {1,2}; each step uses an IMH latent move, the frozen inactive omega, and the tree step skipping the inactive row.

1. **Confirmed: a row entering at a fractional shape and moved by the step from its stale omega is exact under neither contract.**
   - The unextended target defines omega_j only while row j is active. The stale omega_j was drawn under an earlier theta, so after entry (theta, omega) is not distributed as pi, which Thm (sec:setters)(a) needs.
   - One IMH step leaves PG(b_j, psi_j(theta)) invariant but does not restore it from that start. So a contract-C caller fares no better here than a contract-I one.
   - In the toy, an activity update with omega integrated out puts the stationary (theta, active) marginal 0.0058 from the target in total variation.
   - Through setData it is worse. LogisticResponse::setData (and nbinom's) calls coldStart, which resets every row's omega to w/4 (nbinom: (y + r)/4), and it clears the mask.
   - So on main the step after setData starts from a cold value at every fractional row, persisting rows included. Contract C is not exact there either. That conflicts with the note's "from the value the row had" and with dec-B607's "served by the step from their installed values", unless the build makes setData keep omega at rows that persist.

2. **The pseudo-prior extension: confirmed, with two corrections.**
   - Extend by g_j = PG(b_j, 0) at each inactive row. Then:
     - the extended target's marginal is the target, since g_j is a proper density;
     - the tree conditional is unchanged, because g_j has no theta and the empty-leaf rule counts rows regardless of the mask (setActiveRows help);
     - the frozen omega is invariant, since the identity kernel leaves every distribution invariant.
   - The toy confirms it: the extended rule's residual is 3e-17 and its TV is 2e-16, and the latent and tree steps each leave the extended target invariant.
   - **The logistic factor is right.** 2^-b exp(kappa psi - omega psi^2/2), integrated against PG(b, 0), gives back e^{w y psi}/(1 + e^psi)^w, using E e^{-omega psi^2/2} = cosh^-b(psi/2).
   - **Correction: kappa.** With the engine's y in [0, 1] and weight w, the factor needs b = w and kappa = w(y - 1/2). Report5's "y_j - b_j/2" holds only if y_j means the count.
   - **Binomial:** the same factor, times the caller's C(n, k)^w, which the engine never carries.
   - **nbinom:** b = y + r, kappa = (y - r)/2, psi = eta - log r, times Gamma(y + r)/(Gamma(r) y!).
   - **Correction: nbinom with r estimated breaks the cheap route.** g_j then depends on r, so the r conditional gains a factor from each inactive row (the PG(y_j + r, 0) density at its frozen omega), which the engine's r step ignores. The cheap route needs a fixed r. With whole-number y and r the shapes are whole numbers and fresh draws already make it exact.

3. **Confirmed: a caller can implement it on main.**
   - getLatents returns the stale value at an inactive row, as its help says.
   - getFitsWithoutOffset reports f(x_i) at inactive rows; psi_j is that plus the caller's offset.
   - setActiveRows draws nothing.
   - The engine skips inactive rows in refreshLatents, and setState keeps a row's omega across masks.
   - The setActiveRows help's "the latents integrated out ... need not call setResponse" pattern is the inexact one at fractional shapes, so that sentence needs changing.
   - A start is not exact: a row inactive from creation holds w/4, not a draw from g_j. That is a start effect, not a loss of invariance.

4. **Limits: confirmed, with setData wider than stated.**
   - **A weight drawn from 0** needs the ratio PG(omega | w', 0)/PG(omega | w, 0), a density evaluation. Confirmed.
   - **setData** breaks the extension at every inactive row on every call, not only at rows it adds, because coldStart overwrites them all and the mask is cleared.
   - **A fixed-weight toggle** works only if the weights plan's weight 0 skips the draw as the mask does. On main a zero count is refused for logistic, and refreshLatents draws lround(w) PG(1) terms, so fractional shapes exist only in the weights plan.

Reading and seconds-long runs only; nothing written outside the referee's folder.