I confirm the critic's main claim and refute one of its sub-claims. A cheap exact route exists for a contract-C caller when a row is reactivated by setActiveRows or toggled between weight 0 and a fixed weight, but not for rows setData adds.

1. **Confirmed.** Theorem (sec:setters)(a) assumes a fixed set of rows, all carrying omega in the target. When the active set A is part of phi, the joint target only defines omega on A(phi). So a move that brings row j in cannot condition on omega_j: on the unextended target it is necessarily a contract-I move for row j, whatever the caller does with the other rows.
2. **At a fractional weight that is case (c).** The marked row is served by the step from its stale value and is not exact until the step forgets it. The dbartsSampler help's own pattern, "a caller that draws the mask from the fit with the latents integrated out", is exactly this. At whole-number weights the fresh draw makes it exact, by (b).
3. **Refuted: "no extension of the target to inactive rows survives the tree step."** That holds only for a theta-dependent extension. Extend the target with a reference density g_j(omega_j | phi) for each inactive row that does not depend on theta, as in Carlin and Chib's (1995) pseudo-priors. Then:
   - the tree conditional is unchanged;
   - the frozen inactive omega (the engine skips those rows, per model.hpp) is invariant for g_j;
   - the marginal of the extended target is the target.
4. **The cheap rule: take g_j = PG(b_j, 0), with b_j the row's shape when active.** For a logistic row the active factor divided by the inactive one is then 2^(-b_j) exp(kappa_j psi_j - omega_j psi_j^2 / 2), with kappa_j = y_j - b_j/2. That is the same augmented-likelihood term a contract-C caller already uses.
   - A caller exact under the extension updates activity from its conditional given every row's omega, the stale ones included, using that factor for each row whose activity changes.
   - The tools exist on main. getLatents returns the stale value at an inactive row (the help says so), and fits are kept at inactive rows.
   - Engine cost: none, as long as inactive omega are never touched. That is true today, and setState carries them.
5. **Where the rule applies.** It covers setActiveRows reactivation, and a weight toggled between 0 and a fixed w_j (pseudo-prior PG(w_j, 0)).
6. **Where it does not.** It does not cover a weight that moves from 0 to a value the caller draws, because the ratio then needs PG(omega | w', 0) / PG(omega | w, 0), a density evaluation. Nor does it cover a row setData adds, whose omega starts as creation starts it, not as a draw from g_j. Exact draws of PG(b, 0) at a fractional b do not exist either.
7. **Recommended text.**
   - **Theorem (a):** add "for a fixed set of rows; an activity change at a fractional shape is covered only on the pseudo-prior extension".
   - **Critic's help sentence:** keep it for added rows and for weights that move from 0. Replace it for setActiveRows and fixed-weight toggles with the pseudo-prior route and a two-line example: a membership step using exp(kappa psi - omega psi^2 / 2) from getLatents and getFitsWithoutOffset.
8. **Checks.** Reading only: the note's sec:setters, dec-B607, model.hpp on main, the dbartsSampler help on main, and amendment 9 of the critic's report. No runs, no files written.