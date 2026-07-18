# r-ingestion-cleanups

agent: Sonnet
rng: neutral (error paths and dedup only)
budget: ~200 lines

## Goal

Close the R-surface findings of the 2026-07-17 data review that the
remediation arc (data-review-remediation.md, mutator validation) does
not cover.

## Findings

1. Classification-family routing duplicated with drifted wording
   across dbarts.R:294-333, xbart.R:78-113, rbart.R:332-358 (the
   2-vs-3-level check, auto->probit resolution, conflict errors; only
   classifyResponse/announceAutoFamily are shared). Fix: one
   resolveClassificationFamily(data, family, caller) helper; align
   the messages.
2. Sparse test-data asymmetry: dbartsData validity accepts x.test as
   matrix or dbartsMixedMatrix only (A_class.R:482-491), so
   validateXTest silently densifies a bare dgCMatrix test set
   (data.R:112-128) while mixed-container test data stays sparse; no
   test covers the coercion. RESOLVED (VD 2026-07-18): accept
   dgCMatrix test input SYMMETRICALLY - a bare sparse x.test flows
   like the mixed container's, never silently densified. Add tests
   for the sparse-test path (train sparse + test sparse; train dense
   + test sparse behavior defined and tested, erroring informatively
   if the train store cannot serve it).
3. The three copy-pasted NA-policy guards (dbarts.R:962-966,
   997-1001, 1010-1014): shared checkMissingPolicy() helper;
   cosmetic, ride along.

## Constraints

- Suite + all equivalence anchors bitwise; air format --check . and
  lintr clean; error messages keep the package's voice (sentence
  case, no periods per existing style - check neighbors).
- Verify each duplicated site's behavior is truly identical before
  unifying; any intentional divergence gets kept and documented, not
  averaged away.
