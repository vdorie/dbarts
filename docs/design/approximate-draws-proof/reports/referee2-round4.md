# Round 4 confirmation: the authors' round-3 answers

**Verdict: all six points ANSWERED.** Two items have small gaps. Neither affects exactness, and both can be fixed without another round.
- **Item 3:** the S1 threshold rule leaves the range above |psi| = 150 undefined.
- **Item 4:** one sentence in Appendix A claims a tighter bound than the charge it describes.

I concede the cost point: my round-3 reading that the route wins at 1e3 and 1e5 at |psi| = 150 was wrong. I scaled the design's measured microseconds by the old node ratios. Those microseconds were measured at small K, so the scaling hid the K+1 terms each q node costs at |psi| = 150 (K = 223 to 560). With absolute counts the route loses there, as the authors say.

**1. ref1 R1, the cache keyed on the point. ANSWERED.**
- An entry stores (x, b, psi, enclosure). It is used only when x equals omega_i bit for bit and b is unchanged; otherwise the step evaluates fresh. Both cache writes in Algorithm 2 store the tuple.
- The entry is a valid enclosure of log f at (x, b, psi), whatever the history. Eq. (4) moves it exactly to the current psi. So every read gives a valid enclosure, and Theorem 6.1's "decisions taken are exact" holds.
- That covers the four writers named: setData, setState, a fresh draw by the sum, and an outer step that writes the latents.
- The refinement line reads only the enclosure set at initialization. Initialization is the one read site, and it is keyed.
- The same rule is stated in Sec. 1, the text after Algorithm 2, Sec. 15 S2 and Sec. 16 items 2 and 7.

**2. My point 2, "three fresh evaluations". ANSWERED.**
- The note says three only after a RowStep that left a valid entry.
- It says four otherwise: after a fresh draw (a marked row, the first sweep, a row the rule moved from the sum, or nbinom after r moved) or any other write. It gives about 38 percent of steps in the design's estimated-r cell.
- It notes that the design's cache had the same gap.

**3. Cost (my points 3 and 4, ref1 R2). ANSWERED. The table holds; the pass rule is sound and exact, with two gaps.**

Cells I checked, all on the a-priori L10 setting. The step's time is the note's own pricing (116 ns per f node or cell, 12.6 and 3.1 ns per q node and cell term) applied to the counts:

| Cell (b, psi, z) | How checked | f terms | q terms | Density part of a step | Today's sum |
|---|---|---|---|---|---|
| (16, 150, 0) | my own round-3 run, before cost3.py existed (f 102 nodes + 80 cells; q 118 + 80 at K+1 = 97) | 182 | 19206 | 357.7 us | 2.2 us |
| (64, 8, 0) | reran cost3.py (ONLY=L10) | 181 | 985 | 38.2 us | 8.7 us |
| (1000, 50, 0) | reran cost3.py (ONLY=L10) | 238 | 16348 | 346.8 us | 136 us |

- All three match the table exactly: 3.86e4 terms and 358 us; 38.2 us; 347 us.
- The crossings read from the table hold:
  - |psi| <= 8: about 220 to 260 at psi = 0 and about 300 to 400 at psi = 8;
  - |psi| = 50: the route loses at 1e3 (347 against 136 us) and wins at 1e4 (951 against 1360 us);
  - |psi| = 150: it loses at 1e3, 1e4 and 1e5 (1.66 against 0.136 ms; 5.1 against 1.36 ms; 23.5 against 13.6 ms).
- The losses at |psi| = 150 are robust to the pricing. Plain-double node times without radius bookkeeping are lower estimates, and pricing f cells as nodes is an upper one. The crossings at |psi| <= 8 are uncertain by perhaps 2x, which is what S1's timing is for.
- Charging final-h counts rather than the walk is the build Cor. 5.10(b) describes, and it is stated.

The pass rule:
- **Exact.** Under Thm 10.1, chi_i = 1 iff b_i <= b*(|psi_i|) is a function of (b_i, psi_i), not omega_i. Each kernel leaves PG(b_i, psi_i) invariant.
- **Sound as a decision procedure.** b* is the smallest whole b with c_R <= c_S, taking the larger bracketing value. The gate is b* <= 1e3 at |psi| <= 8; failing it means no whole-shape case, and S1 goes back to the maintainer before S2.
- **Gap (a), minor.** "The larger b* of the bracketing grid points" leaves |psi_i| > 150 undefined. State it: b* = b*(150), which is infinity in the scripts, or infinity outright.
- **Gap (b), minor.** c_R should be the expected per-row cost, including the share of four-evaluation steps. That share is 38 percent in the nbinom cell, plus rows that the psi-dependent rule moves between the sum and the step.
- **Gap (b), continued.** S1's timing grid (b in {1, 16, 64, 1e3, 1e4, 1e5}) has no fractional shape, yet the 70 us budget for real weights needs one. Add b = 0.2 and 2.5.

**4. My point 1 and ref1 R3, App. A's l_f bound. ANSWERED, with one inconsistent sentence.**
- These are all correct:
  - R_W is charged explicitly;
  - the count is gamma^(2)_12;
  - log cos v is formed as log sin(pi/2 - v), with pi/2 in double-double;
  - (A2') names sin, and eq. (8) and Q1 list sine;
  - the saddle relation b tan(v)/(4v) = x, and |W0| <= |t0|/2 when W0 < 0;
  - the leftover error from the log sin form, about pi^2 x u_2 <= 2 u_2 x |t0|.
- **The inconsistent sentence.** On W0 >= 0 the charge is b u_2 |W0| max{1/2, tan(v)/(2v)}. Continued to v = i w0, the second term is tanh(w0)/(2w0) < 1/2, so the charge is b u_2 W0 / 2 = b u_2 w0^2 / 2. The text says this is at most b u_2 w0 / 2, which is false for w0 > 1.
- **Numbers at (1e5, 150).** The charge is about 5.5e-23, and the main term (gamma^(2)_12 + 2 u_L2)(1.49e7) is about 5.9e-23. Together that is about 1.1e-22, a hair above the stated "below 1e-22". This is negligible against the 1e-12 floor.
- **Fix.** Charge b u_2 |W0| tanh(w0)/(2w0) on W0 >= 0, which is at most b u_2 w0 / 2. Or say "about 1e-22".

**5. My point 5, the update's operation count. ANSWERED.** Sec. 1, item 7 and Q1 all say two double-double exponentials and log1p's, about ten double-double operations, a few hundred double operations. Q1 adds q's 2(K+1) double-double log1p's.

**6. What round 3 broke. Nothing found.**
- Lemma 5.9, checked line by line, holds:
  - the expansion of abs(z^2 + tau'^2)^2;
  - g(A) = 2(A - 2u^2)/A^2 decreasing on A >= 4u^2, so kappa = g(min A_j);
  - the closed forms Prod A_j(u)/A_j(rho) = cosh sqrt(W_u) / cosh sqrt(W_rho), T_1 = tanh(r)/(2r) and T_2 = (tanh r - r sech^2 r)/(4r^3), with T_1(0) = 1/2 and T_2(0) = 1/6.
- Cor. 5.10 holds:
  - (a): the J_0 terms do not increase for Y >= A_j - 2rho^2; the delta condition; the threshold for (rho + y) e^{-delta y^2} to decrease; the sum-to-integral bound; the integral bound (rho/y1 + 1)/(2 delta);
  - (b): the h formula.
- agg_check.out: the bound held at every point (max lhs - rhs = -0.0038).
- cacheupd_r3.out: every moved enclosure overlaps a fresh one. At (1e5, 150, 3) with dpsi = 5 the midpoints differ by 9.1e-13, inside the combined radius of 1.4e-12.

Not checked: cost cells other than those three; termbench.c and pgbench.c themselves; ref1's R4 and R5 beyond reading the text. I ran only seconds-long probes this round, wrote and edited no files, and left nothing running in the background.