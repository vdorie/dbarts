# Approximate draws: Polya-Gamma latents whose exact cost grows with their shape

Status: PLANNED 2026-11-07, on the approximate-draws branch, not landed; slice S1 (the certified density and one
row's step, wired into nothing) built
([approximate-draws: Polya-Gamma latents whose exact cost grows with their shape](../plans/approximate-draws.md#approximate-draws-polya-gamma-latents-whose-exact-cost-grows-with-their-shape)).
The proof and the S1 note beside it: [The accept-or-keep Polya-Gamma step: proof](approximate-draws-proof/README.md#the-accept-or-keep-polya-gamma-step-proof).

The class: S1 neutral (no draw of any existing call moves); the exact route shifting where it is selected; real
shapes a new capability; the approximate settings posterior-changing where a user selects them.

## The question

Each row of a logistic, binomial, logistic-hazard or negative binomial fit carries a Polya-Gamma precision
omega ~ PG(b, psi), redrawn every sweep. Today it is drawn exactly as the sum of b draws of PG(1, psi), so its cost
grows with the row's shape b (trials times weight, or count plus dispersion), and a real-valued b has no draw at
all. How is omega drawn at a cost that does not grow with b, exactly where that is possible, and approximately
only when a user asks?

## The exact route

At the start of each sweep, before the trees, a row takes an accept-or-keep step: it proposes omega' from q, a sum
of K + 1 gammas, and keeps omega' with probability min(1, r), r = f(omega') q(omega) / (f(omega) q(omega')),
otherwise keeps omega. f is PG(b, psi) and q the density of what the sampler actually draws: K gammas Ga(b, r_k)
at the target's first K rates, r_k = 2 pi^2 (k - 1/2)^2 + psi^2 / 2 formed as doubles (the first rounded down so
that r_1 < d_1), and one Ga(alpha, beta) matching the mean and variance of the remainder. The step leaves
PG(b, psi) invariant for every b > 0, psi and K >= 1, and rounding r_1 down makes it uniformly ergodic; the sweep
then leaves the joint posterior invariant. That is the proof's result ("The row kernel", "The sweep"); it is the
package's own construction, so it was proven by independent agents before the build.

K comes from a rule on D, the proposal's distance from f in closed form: the smallest K with b D(psi, K) below a
fixed tolerance ([`ruleK`](../../src/bartcore/pgdensity.hpp)). The proposal is
[`Proposal`](../../src/bartcore/pgdensity.hpp).

## The density, certified

The step needs log f and log q at two points each, and a decision taken on them must be the decision exact
arithmetic would take. Each density is evaluated as an enclosure: a midpoint and a radius, the true log density
within the radius, every rounding charged. One rule serves f and q at every shape and x > 0: the inverse Laplace
transform on the parabolic contour through the saddle point, by the trapezoid rule with an a priori bound (edge
lines on either side bound the discretization error; each dropped node's majorant and a monotone tail bound the
truncation). The node exponents are formed in a cancellation-free form relative to the saddle, and the lead term
x t_0 + log L(t_0) in double-double, so the radius stays near 1e-14 of the log density's scale at large shapes.

[`TargetDensity`](../../src/bartcore/pgdensity.hpp) is f, [`GammaSumDensity`](../../src/bartcore/pgdensity.hpp)
is q (any sum of gammas); [`logDensity`](../../src/bartcore/pgdensity.hpp) evaluates either at a refinement level.
Double-double, with five functions whose relative errors are stated and tested against 200-digit references, is
[ddouble.hpp](../../src/bartcore/ddouble.hpp).

The certificate comes in three levels fixed in the code ([`kLevelsF`](../../src/bartcore/pgdensity.hpp),
[`kLevelsQ`](../../src/bartcore/pgdensity.hpp)): a cheap one (one edge line a side, at most 40 cells a line), a
middle one and the full default. Each is a complete bound; the step evaluates the cheap one and refines only when
the decision is not yet taken, and intersects the enclosures it has.

## The decision and close calls

The uniform is held as an interval of width 2^-53 from two 32-bit words of the chain's generator. The step accepts
when the enclosure of log r lies above the interval's upper end (or above 0), keeps when it lies below the lower
end, and otherwise refines. A decision taken is the one the exact uniform and exact r give. After the last level a
step still undecided is a close call: the row keeps its value, and the call is counted (S2 reports the count and
warns). The chance of a close call at a step is at most the floor, twice the enclosure's half-width (the four
densities' radii and the sum's rounding) plus the thresholds' slack and the uniform's grain
([`stepFloor`](../../src/bartcore/pgdensity.hpp)); S1 measures it at the pins at most 9.6e-12, against the gate's
2.3e-11. [`stepRow`](../../src/bartcore/pgdensity.hpp) is the one-row step.

## Where the step costs less than the sum

The step's cost does not grow with b, the sum's does. A whole-number row takes the step only above a threshold
b*(|psi|) at which the step measures twice the sum's speed, never below 64, and the sum otherwise; a row with a
fractional shape always takes the step. S1 times both on one core; the threshold's dependence on |psi| goes to the
maintainer with those timings before S2 uses it. Above |psi| 150 the sum is always used.

## What each slice builds

- S1: the certified densities, the proposal and the one-row step, unused by the sampler; their tests in
  tests/cpp (pins against reference tables generated in Arb, the double-double gate, the K table, the floor, the
  reach into both tails, the step's distribution); the timing harness, benchmarks/kernels/pg_step.cpp.
- S2: the step in the sampler for logistic, binomial, logistic-hazard and negative binomial rows, the per-row rule,
  the marks a setter leaves, the close-call count.
- S3: real weights and shapes under the step, every setter taken at fractional shapes (dec-B607).
- S4: the series with a moment-matched remainder as an approximate method, chosen by a tolerance.
- S5: the method by name at the doors users reach.

## The proof

The proof note, its scripts and the referee reports are in [The accept-or-keep Polya-Gamma step: proof](approximate-draws-proof/README.md#the-accept-or-keep-polya-gamma-step-proof):
the note gives the step's exactness, the density's a priori bound and its rounding in ball arithmetic. The S1 note
there, s1-evaluator.tex, states and proves what the engine's evaluator adds to the note's: binary64 ball rules,
the lean node form for f at large shapes and the second-order form for q, a smaller J_0 in the aggregated
majorant, the monotone tail's constants, log sinh past a double's range, the dyadic far tail for q, the edge lines'
distance cap and the uniform as an interval, with the cheaper certificate's settings.
