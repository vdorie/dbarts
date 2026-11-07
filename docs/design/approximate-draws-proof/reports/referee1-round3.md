# Referee 1, third round: the revised accept-or-keep Polya-Gamma note (proof.pdf)

**Verdict: minor revision.** The math in the revision holds: the cache update, the double-double lead and its bound, the new Prop. 6.2 constant and the Prop. 6.4 proof. Two things need fixing before the note can stand:
- **A one-line spec fix.** The cache is keyed on (b, psi) only, so it goes stale when omega_i is replaced by anything other than the step (R1).
- **A restatement of the cost conclusion.** Calling the cost open is right, but the note understates what is at stake. The count ratios leave out the K+1 terms per q node. On the authors' own counts, at |psi| = 150 the route is already near or above today's exact sum at b = 1e3 to 1e5 (R2).

The rest are minor (R3 to R5).

## 1. Round-2 findings

**N1 (cache; four evaluations a step): resolved for the count; the repricing is partly resolved (R2).**
- **The identity is correct.** Eq. (4) is the tilt formula at psi' divided by that at psi. lambda is even, so a sign change of psi is fine.
- **I reran it at new points** with cacheupd.py, update at 106 bits. Every update contained the fresh evaluation, with update radius 8e-28 to 1e-31:
  - (3e4, psi = -60, z = 5) with dpsi = 0.7, -2.3 and +121 (a sign crossing);
  - (0.5, 10, z = 8) with dpsi = 0.3 and -10.5.
- **The count is right.** A step makes three fresh evaluations, f(y), q(y) and q(x), plus the update. q(x) must be fresh, because the proposal is re-matched at psi.
- **New gap.** The cache is invalidated only on a change of b (R1).

**N2 (floor assumes a 106-bit lead): resolved.**
- The lead is now specified in double-double (eq. (8)) and given a rounding bound (App. A). The bound uses u2 = 16u^2 = 2^-102, which covers the Joldes-Muller-Popescu bounds for accurate addition (about 3u^2), multiplication (5u^2) and division (15u^2 + 56u^3).
- I checked the arithmetic of the bound:
  - f: 1.5e7 x (gamma_10^(2) + 2 u_L2), about 5e-23;
  - q: 1.1e5 x (gamma_564^(2) + u_L2), about 1e-23.
  - Both are far below the node radii, so the 100-bit stand-in is adequate whatever their exact order.
- The plain-double form is rejected, with a certified 2e-9 that fits my measured 9.5e-11 as a worst case.
- The engine's floor is now an S1 gate. One justification is wrong but the bound survives (R3).

**N3 (branch cuts): resolved.** I checked the new rule:
- the upper-half-plane facts hold, since Im W = 2 rho y >= 0 and 2 tau_j^2 + t0 > 0;
- the principal sqrt and log agree with their limits from above on the cut, so the intersection rule is valid;
- the fallback to the direct form for the log1p of eq. (7) is sound.

**N4 (Prop. 6.4): resolved.**
- The two-state proof is correct.
- I checked Toy 6's two-state numbers by hand, with target (2/3, 1/3), uniform proposal and d = 1/100:
  - the flows are 0.1633 and 1/6, so the residual is 1/300;
  - the stationary distribution has pi_0 = 1/1.49, so the total variation is 0.00447.
- The capped-ladder paragraph is right.

**N5 (Prop. 6.2's constant): resolved.**
- p_j^{3/2} is right: h_j is of order 1/p_j, and under e^{-2xy^2} decay the truncation point is of order sqrt(p_j).
- The recursive-sum rounding is gamma_N times the sum of the moduli, and that sum converges.
- The polynomial-factor remark closes the cost argument.

**N6 (Thm 11.2(c) too broad): resolved.** The narrower wording is in the theorem, T7, Algorithm 1 and Q2.

**N7 (Cor. 11.4's construction): resolved.**
- The quantile rotation is a contract-I step in every model with a continuous conditional.
- The countable dense set plus continuity of the Laplace transform gives Lemma 11.3's hypothesis.
- psi_0^2/2 sweeps an open interval even when psi_0's range contains 0.

**N8 (pin maxima; shared q reference): resolved.**
- "At the pins" is added.
- The q pins now have two references independent of the contour:
  - the vertical line, where the total shape is large;
  - Moschopoulos's series, where the largest rate times x is small.
- In pins_q_B2.out they agree with the contour references to 2e-35 to 1e-17, with act/rad at most 0.53.
- The vertical line's failure at b = 0.2, |psi| = 8 is shown in pins_q_B.out (refdis up to 1). It is the reference's failure, as the note says.

## 2. New material

- **Cache (eq. (4)): correct, but incompletely keyed (R1).**
  - Its double-double radius is below 1e-23, as stated.
  - Storing the cached value in double-double is needed and is specified.
  - The authors' check script adds the update to the cached value at 53 bits, which is a script artifact (R5).
- **Lead (eq. (8)) and App. A: correct.** One justification fails on the W0 < 0 branch, though the bound survives (R3).
- **Floor quoted for the double-double lead: acceptable.**
  - 6.15e-12 + 2^(1-beta) is the right arithmetic: 2(2 x 1.41e-12 + 2 x 1.28e-13).
  - The 100-bit stand-in shows the same maxima over 120 f pins and 30 q pins.
  - Ordering Arb at 100 bits against App. A's bounds is plausible but unproven. It does not matter at 1e-22 against 1e-12.
- **S1 gate at 2.3e-11: sound as a gate.**
  - It is 2^-32/10, which leaves about 3.7 times Arb's figure for the engine's slack.
  - It is still at least 20(K+1) times below the error the step's own uniforms and gammas carry.
  - It sits awkwardly with Q1's user-facing figure (R4).
- **Prop. 6.2's C p_j^{3/2} 2^{-p_j}: correct** (N5 above).
- **Cost repricing (Sec. 16 item 8): the conclusion "open until S1 times it" is right as far as it goes, but understated.** The counts' method has three defects (R2).

## 3. Findings

**R1. Gap (correctness of the specification). The cache must also be keyed on omega_i.**
- **What the note says.** The cache holds log f at the current point with "the (b, psi) it was computed at", is "moved when psi_i changes and recomputed when b_i changes" (S2, Sec. 15; Sec. 16 item 7). Algorithm 2 initializes I_f(x) from the cache unless there is none or b changed.
- **Where it goes stale.** omega_i can be replaced without b_i or psi_i changing, and the stale value then gives a wrong log r:
  - **setData** marks every row, and marked integer-shape rows at b > 64 are drawn fresh. With only y changed in a logistic or binomial fit, (b_i, psi_i) is unchanged, since PG(b, eta) does not depend on y.
  - **setState** installs omega.
  - **Contract-C outer steps** write the latents.
- **Consequence.** The step is then not invariant, silently.
- **The old design** had the same implicit assumption: its log f0(omega) cache.
- **Fix.** Store the exact double x with the cache entry and use the entry only if it matches omega_i. Or invalidate the entry on every write to omega_i outside the step.
- **Evidence.** Reading of Sec. 1, Algorithm 2, Sec. 15 S2, Sec. 16 item 7 and Thm 11.2. A grep for invalidation and recomputation finds only the b clause.

**R2. Gap (cost; the conclusion is understated). The ratios leave out K, the absolute counts already compete with today's sum, and the counts cannot be reproduced.**
- **(a) K+1 is hidden.**
  - The table counts a q node or cell as one unit, "as the design's q nodes do". So the ratio is a time ratio only where K matches the sizing's fits.
  - The rule's K at the table points: 1 and 37 at b = 0.2; 1 and 96 at 16; 3 and 223 at 1e3; 6 and 560 at 1e5, each for |psi| = 0 and then 150 (kd.K_rule).
  - Term evaluations per step, f + 2q(K+1), from the authors' counts:
    - (16, 150): about 1.5e4;
    - (1e3, 150): about 1.0e5;
    - (1e5, 150): about 2.1e6.
  - Today's exact sum costs 160 ns x b (design table): 160 us at 1e3 and 16 ms at 1e5.
  - At 5 to 20 ns per term (my assumption, not measured), the route's density part costs:
    - at (1e3, 150), 0.5 to 2 ms: slower than the sum even at the low end;
    - at (1e5, 150), 10 to 42 ms: comparable to the sum.
  - So the stakes are more than "the b > 64 threshold may move". At large |psi| the certified route as specified may not beat the sum it replaces at any shape. dec-B590's rule may need K or |psi| in it, not b alone.
- **(b) The counts cannot be reproduced and are tuning-sensitive in both directions.**
  - "The best of three cell settings" are LEAN settings that are recorded nowhere, and there is no output file.
  - My own runs (chk3/lean.py):
    - **(16, 0, z = 0), default settings:** 8410 node-equivalents for f against the table's 110, a factor of 76.
    - **(16, 0, z = 0), my lean setting** (one edge line each side, 10 to 40 cells, no splitting): 212. Its discretization bound is 2e-23 relative, about a billionth of the rounding radius.
    - **(1e5, 0, z = 0), default settings:** f costs 125607 node-equivalents. The discretization bound is 4.7e-83, so most of that work buys nothing.
    - **(1e5, 0, z = 0), four lean settings of mine:** every one gave an infinite discretization bound.
  - costcount.py also prints the radius from the default run, not the lean run's. So the table does not show that its own counts reached the floor.
  - Also, the counts are at z = 0 only.
- **(c) What S1's timing settles.** S1 timing the certificate as specified will reproduce roughly these counts. What is actually open is a cheaper certificate:
  - the combined edge-line majorant;
  - a monotone tail in place of the node-by-node truncation, which accounts for about 900 of the roughly 1000 node evaluations for f at (1e5, 0) in my run;
  - a coarser h, since the floor-level h overshoots the rounding radius by about 1e6.
  - That is design and proof work. It should be specified before S1, or made part of S1, and then timed.
- **(d) Wording.** response2.md and final2.md label the ratio "per-row cost". The note, correctly, says "the density part"; the ratio leaves out the K+1 gammas and the saddle solves. The summary's "3 to 16" should read 2.5 to 16.
- **Fix.**
  - Record the LEAN settings and the output, and print the lean run's own radius.
  - Add an absolute column (term evaluations) and today's sum at the same b.
  - State that at large |psi| the route may lose to the sum.
  - Make the cheaper certificate an S1 item, ahead of timing.

**R3. Minor. App. A's bound for the radius of l_f gives a wrong reason on the W0 < 0 branch.**
- **The problem.** The bound is justified by "the rounding of W0 being carried through lambda with |lambda'| <= 1". On the log cos branch, d/dv log cos v = -tan v is unbounded.
- **The bound still holds.** At the saddle, b tan v = 4vx. The term from the rounding of W0 and v is then about 4x|W0| u2 <= 2x|t0| u2, since |t0| >= 2|W0|. That is covered by the |x t0| term.
- **Two edits needed.**
  - Use the saddle relation to justify the bound. It holds only near the saddle, so either state that or charge b v tan v u2 explicitly.
  - Make (A2') say that cos has a relative error bound near pi/2. Or form log cos v as log sin(pi/2 - v), with pi/2 held in double-double.
- **Evidence.** Reading of App. A and eq. (8).

**R4. Minor. The S1 gate against Q1's user figure.**
- **The mismatch.** Q1 says "fewer than one in a thousand fits of 1e4 rows by 1e4 sweeps". At the gate's limit of 2.3e-11, one fit's 1e8 row steps give about 2.3e-3, which is one fit in about 430.
- **The fix.** Q1 should tie the user figure to the engine's measured floor, or quote the figure at the gate.
- **A gap.** The note does not say what follows a failed gate: option B', or a tighter evaluator.

**R5. Minor (script).** cacheupd.py forms lf = r0['logp'] + upd at 53 bits. Its "updated rad" therefore includes a 53-bit rounding of the sum: 1.06e-13 grows to 2.19e-13 at (3e4, -60), dpsi = -2.3. The engine spec stores and adds in double-double, so the claim stands, but the script does not check the double-double addition. Fix: add at 106 bits.

## 4. Not checked

- The 150 pins at 100 bits individually: I read the maxima only.
- The authors' own costcount table entries: their LEAN settings are not recorded.
- Per-term timings: my 5 to 20 ns is an assumption.
- Arb's correctness, which I assumed.
- Sections not changed since round 2.

My script: the referee's own scripts (not kept) It runs costcount.py's count under a given cell setting and also reports the radius that setting itself gives. Reruns used the authors' scripts/cacheupd.py and costcount.py and kd.K_rule; I ran nothing else. One core, about 45 minutes; no background jobs.