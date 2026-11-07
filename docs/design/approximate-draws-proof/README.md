# The accept-or-keep Polya-Gamma step: proof

"The accept-or-keep Polya-Gamma step: where it keeps the posterior exact, and where it does not" (2026-11-06),
checked by two independent referees. Its source is `texparts/p1.tex` to `p4.tex`; `build.py` joins them, fills
the numbers from `scripts/fill.txt` and the tables it names, and runs pdflatex twice. No PDF is checked in.

```sh
python3 build.py OUTDIR    # writes OUTDIR/proof.tex and OUTDIR/proof.pdf; needs pdflatex on PATH
```

## The engine's evaluator

`s1-evaluator.tex` is a standalone note: what the engine's certified evaluator (src/bartcore/pgdensity.hpp) adds
to the proof's specification, each departure proven, and the cheaper certificate's settings.

```sh
cd OUTDIR && pdflatex s1-evaluator.tex && pdflatex s1-evaluator.tex
```

## Scripts

`scripts/` holds the checks the proof's appendix "Numerical verification" lists, each with the `.out` file it
printed. They need Python 3 and the packages in `requirements.txt` (mpmath and python-flint, for Arb):

```sh
python3 -m venv VENV && VENV/bin/pip install -r requirements.txt
cd scripts && ../VENV/bin/python exact_toys.py
```

Each script's docstring gives its usage; the scripts import one another from `scripts/`, so run them from there.

`scripts/gen_pgdensity_refs.py` writes the engine tests' reference tables, tests/cpp/pgdensity_refs.inc (about
seven minutes on two processes):

```sh
cd scripts && ../VENV/bin/python -B gen_pgdensity_refs.py ../../../../tests/cpp/pgdensity_refs.inc 2
```

`scripts/bench/` times today's exact sum and the kernels' per-term costs. From the package root, with R's headers
and library:

```sh
B=docs/design/approximate-draws-proof/scripts/bench; RH=$(R RHOME)
cc -O2 -Isrc/include -Isrc -I$RH/include -o pgbench $B/pgbench.c $B/stubs.c src/external/random.c \
  src/external/randomBase.c src/external/randomNorm.c -L$RH/lib -lR -lm
cc -O2 -o termbench $B/termbench.c -lm
```

## Referee reports

- `reports/referee1-round3.md`: the first referee's last round on the note.
- `reports/referee2-round4.md`: the second referee's confirmation of every point of that round.
- `reports/referee2-reentry.md`: rows that re-enter the fit, the source of the pseudo-prior route.
- `reports/referee1-reentry.md`, with its exact toy `reports/referee1-reentry-pseudo.R`: the first referee's
  check of it.
