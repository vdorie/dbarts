#include "common.hpp"

// ---------------------------------------------------------------------------
// Property fuzzer for the transactional mutation surface. A seeded ext_rng
// drives long random op sequences over mixed sampler configurations; after
// every op the same invariants are re-checked, so rollback and re-routing
// bugs surface with a replayable seed and op trace.

static double fuzzUnif(ext_rng* r) {
  return ext_rng_simulateContinuousUniform(r);
}

static size_t fuzzInt(ext_rng* r, size_t k) {
  return k <= 1 ? 0
                : static_cast<size_t>(ext_rng_simulateIntegerUniformInRange(
                    r, 0, static_cast<int64_t>(k)));
}

// Owns every buffer handed to the sampler as a borrowed pointer, so lifetimes
// outlive a whole sweep regardless of which op installed them.
struct FuzzArena {
  std::vector<std::unique_ptr<std::vector<double>>> bufs;
  double* keep(std::vector<double>&& v) {
    bufs.push_back(std::make_unique<std::vector<double>>(std::move(v)));
    return bufs.back()->data();
  }
};

// owner[i] <- the bottom node whose index range currently holds observation i,
// read off the live partition alone.
static void fuzzFillOwner(const Tree& t, int32_t i, std::vector<int32_t>& owner) {
  const Node& nd(t.at(i));
  if (nd.isBottom()) {
    for (size_t m = nd.begin; m < nd.end; ++m) owner[t.indices[m]] = i;
    return;
  }
  fuzzFillOwner(t, nd.leftChild, owner);
  fuzzFillOwner(t, nd.leftChild + 1, owner);
}

// The state a rejected or abandoned mutation must leave untouched. The old
// fingerprint carried codes, cuts, sigma and Chain::treeFits(), which is
// forest 0 alone: a rollback that restored mu and left tau routed by the
// proposal passed it green. It is rebuilt here so that omission stops being
// expressible (docs/plans/multiforest-predictor-mutation.md, "The snapshot,
// closed structurally"): the persisted state entire, plus every LIVE structure
// the persisted state does not carry, captured by INDEX LOOPS with no forest
// literal anywhere - a new forest class is covered on arrival. F4
// (testSnapshotCoversEveryFamily) is the gate that each captured family can
// make the comparison false.

// One arena node's structure: the split rule and the index range it owns.
struct FuzzNode {
  int32_t leftChild = invalidNode;
  int32_t variable = invalidVariable;
  std::uint64_t bits = 0;
  size_t begin = 0, end = 0;
  bool operator==(const FuzzNode&) const = default;
};

// One tree's live geometry: its arena, its pooled-mask pool, and the partition
// its index buffer expresses - observation i's owning bottom node, read off
// the live begin/end ranges.
//
// leafOwner, and NOT the raw index buffer: the rollback re-routes from the root
// (Chain::repartitionTrees), and the kernel's partition is not order-stable, so
// a rejected transaction restores each leaf's MEMBERSHIP but generally permutes
// its members. Measured on the shipped rollback, S1 engine and pre-S1 alike:
// every other family below - codes, cuts, the whole persisted state, split
// rules, index RANGES, fits and totalFits - comes back bitwise, and only the
// within-leaf order moves. Recorded here rather than asserted because it is
// pre-existing behavior with a live baseline over it (the predreject
// equivalence scenario runs a second leg after a rejected transaction, so
// whatever the permutation does to a later sweep's reduction order is pinned
// bitwise there).
struct FuzzTree {
  std::vector<FuzzNode> nodes;
  std::vector<std::uint64_t> maskPool;
  std::vector<int32_t> leafOwner;
  bool operator==(const FuzzTree&) const = default;
};

// Everything getState does not carry, per chain and per forest.
struct FuzzGeometry {
  std::vector<std::vector<std::vector<FuzzTree>>> trees;   // chain, forest, tree
  std::vector<std::vector<std::vector<double>>> fits;      // chain, forest
  std::vector<std::vector<std::vector<double>>> totals;    // chain, forest
  std::vector<std::vector<FuzzTree>> varianceTrees;        // chain, tree
  std::vector<std::vector<double>> varianceFactors;        // chain
  std::vector<std::vector<double>> varianceFits;           // chain
  bool operator==(const FuzzGeometry&) const = default;
};

template <typename S>
struct FuzzSnapshot {
  // the store's quantized predictors; the persisted state carries the cut grid
  // but never the codes, and a rollback must put both back
  std::vector<xint_t> codes;
  SamplerStateData state;
  FuzzGeometry geom;
};

static FuzzTree fuzzCaptureTree(const Tree& tree, size_t n) {
  FuzzTree out;
  out.nodes.resize(tree.nodes.size());
  for (size_t i = 0; i < tree.nodes.size(); ++i) {
    const Node& nd(tree.nodes[i]);
    out.nodes[i] = FuzzNode{nd.leftChild, nd.rule.variableIndex, nd.rule.bits,
                            nd.begin, nd.end};
  }
  out.maskPool = tree.maskPool;
  out.leafOwner.assign(n, invalidNode);
  fuzzFillOwner(tree, 0, out.leafOwner);
  return out;
}

template <typename S>
static FuzzSnapshot<S> fuzzCapture(S& s) {
  FuzzSnapshot<S> snap;
  size_t n = s.numObservations();
  snap.codes = s.data().train.codes;
  s.getState(snap.state);
  FuzzGeometry& g(snap.geom);
  g.trees.resize(s.numChains());
  g.fits.resize(s.numChains());
  g.totals.resize(s.numChains());
  g.varianceTrees.resize(s.numChains());
  g.varianceFactors.resize(s.numChains());
  g.varianceFits.resize(s.numChains());
  for (size_t c = 0; c < s.numChains(); ++c) {
    const auto& ch(s.chain(c));
    for (size_t f = 0; f < ch.numForests(); ++f) {
      size_t numTrees = ch.numTreesInForest(f);
      std::vector<FuzzTree> trees(numTrees);
      for (size_t t = 0; t < numTrees; ++t)
        trees[t] = fuzzCaptureTree(ch.treeInForestForTesting(f, t), n);
      g.trees[c].push_back(std::move(trees));
      std::vector<double> slab(n * numTrees);
      ch.forestTreeFits(f, slab.data());
      g.fits[c].push_back(std::move(slab));
      g.totals[c].push_back(ch.totalFitsInForest(f));
    }
    if (!ch.hasVarianceForest()) continue;
    size_t m = ch.numVarianceTrees();
    for (size_t j = 0; j < m; ++j)
      g.varianceTrees[c].push_back(
        fuzzCaptureTree(ch.varianceTreeForTesting(j), n));
    const double* factors = ch.varianceFactorsForTesting();
    g.varianceFactors[c].assign(factors, factors + m * n);
    const double* combined = ch.varianceFits();
    g.varianceFits[c].assign(combined, combined + n);
  }
  return snap;
}

template <typename S>
static bool fuzzSnapshotsEqual(const FuzzSnapshot<S>& a,
                               const FuzzSnapshot<S>& b) {
  // statesAgree covers the per-chain and per-forest half of the persisted
  // state and nothing above it: SamplerStateData's own cutPoints and
  // currentSampleNum have no comparison there (deliberately - the state tests
  // compare states across cut grids), so the sampler-level fields are compared
  // here, where a rollback must restore both.
  return a.codes == b.codes && a.state.cutPoints == b.state.cutPoints &&
         a.state.currentSampleNum == b.state.currentSampleNum &&
         statesAgree(a.state, b.state) && a.geom == b.geom;
}

// Every internal node's children must nest its index range exactly and every
// leaf must be occupied; the running sum then equals n iff the leaves cover it.
static bool fuzzSubtreeCovers(const Tree& t, int32_t i, size_t& leafObs) {
  const Node& nd(t.at(i));
  if (nd.isBottom()) {
    if (nd.numObservations() == 0) return false;
    leafObs += nd.numObservations();
    return true;
  }
  const Node& lo(t.at(nd.leftChild));
  const Node& hi(t.at(nd.leftChild + 1));
  if (lo.begin != nd.begin || hi.end != nd.end || lo.end != hi.begin)
    return false;
  return fuzzSubtreeCovers(t, nd.leftChild, leafObs) &&
         fuzzSubtreeCovers(t, nd.leftChild + 1, leafObs);
}

// Routing agreement: the partition a tree HOLDS must be the partition its split
// rules DESCRIBE. fuzzSubtreeCovers is interval arithmetic over index ranges
// and never consults a rule, so a partition left over from the previous
// predictors satisfies it exactly; this re-descends every observation through
// the rules against the current codes and demands it land in the leaf that owns
// it. That is the invariant a missed re-route violates.
static bool fuzzTreeRoutesCorrectly(const Tree& t, const ColumnStore& data,
                                    size_t n, std::vector<int32_t>& owner) {
  owner.assign(n, invalidNode);
  fuzzFillOwner(t, 0, owner);
  for (size_t i = 0; i < n; ++i)
    if (t.findBottomNodeForObservation(data, i) != owner[i]) return false;
  return true;
}

template <typename S>
static const char* fuzzInvariantViolation(S& s) {
  size_t n = s.numObservations();
  for (size_t c = 0; c < s.numChains(); ++c) {
    const auto& ch(s.chain(c));
    // occupancy, coverage and the totalFits identity, per FOREST: the check
    // read forest 0 against forest 0's total, which is self-consistent and so
    // would not misfire on a BCF - it simply left forest 1 uncovered
    std::vector<double> fits;
    for (size_t f = 0; f < ch.numForests(); ++f) {
      size_t numTrees = ch.numTreesInForest(f);
      const std::vector<double>& total(ch.totalFitsInForest(f));
      fits.resize(n * numTrees);
      ch.forestTreeFits(f, fits.data());
      for (size_t t = 0; t < numTrees; ++t) {
        const Tree& tree(ch.treeInForestForTesting(f, t));
        if (!tree.bottomNodesAreOccupied()) return "empty leaf";
        size_t leafObs = 0;
        if (tree.at(0).begin != 0 || tree.at(0).end != n ||
            !fuzzSubtreeCovers(tree, 0, leafObs) || leafObs != n)
          return "partition does not cover n";
      }
      for (size_t i = 0; i < n; i += 17) {
        double acc = 0.0;
        for (size_t t = 0; t < numTrees; ++t) acc += fits[t * n + i];
        if (std::fabs(acc - total[i]) > 1e-8 * (1.0 + std::fabs(total[i])))
          return "totalFits != tree-order sum";
      }
    }
    if (!(std::isfinite(s.sigma(c)))) return "sigma not finite";
    std::vector<int32_t> owner;
    for (size_t f = 0; f < ch.numForests(); ++f)
      for (size_t t = 0; t < ch.numTreesInForest(f); ++t)
        if (!fuzzTreeRoutesCorrectly(ch.treeInForestForTesting(f, t),
                                     s.data(), n, owner))
          return "mean tree routes an observation to a foreign leaf";
    if (!ch.hasVarianceForest()) continue;
    size_t m = ch.numVarianceTrees();
    for (size_t j = 0; j < m; ++j)
      if (!fuzzTreeRoutesCorrectly(s.varianceTreeForTesting(c, j), s.data(), n,
                                   owner))
        return "variance tree routes an observation to a foreign leaf";
    // s^2(x_i) is maintained incrementally across a sweep and recomputed as the
    // fresh product at its end, so it must equal that product exactly here; a
    // factor is a drawn scale and is strictly positive by construction.
    const double* combined = ch.varianceFits();
    const double* factors = ch.varianceFactorsForTesting();
    for (size_t i = 0; i < n; ++i) {
      double product = 1.0;
      for (size_t j = 0; j < m; ++j) {
        double h = factors[j * n + i];
        if (!(h > 0.0)) return "variance factor not positive";
        product *= h;
      }
      if (std::fabs(combined[i] - product) > 1e-12 * (1.0 + std::fabs(product)))
        return "combinedVariance != per-tree factor product";
    }
  }
  return nullptr;
}

enum FuzzOp {
  OP_SET_PREDICTOR, OP_UPDATE_COLUMNS, OP_PER_OBS, OP_SESSION_ABANDON,
  OP_SET_DATA, OP_SET_CUTS, OP_SET_SIGMA, OP_SET_RESPONSE, OP_SET_WEIGHTS,
  OP_SET_OFFSET, OP_SET_TEST, OP_RUN, OP_STATE, OP_COUNT
};

static const char* const fuzzOpName[OP_COUNT] = {
  "setPredictor", "updateColumns", "perObs", "sessionAbandon", "setData",
  "setCuts", "setSigma", "setResponse", "setWeights", "setOffset", "setTest",
  "run", "state"};

static const int fuzzOpWeight[OP_COUNT] = {5, 5, 4, 3, 2, 3, 2, 2, 2, 2, 2,
                                           3, 2};

static int fuzzPickOp(ext_rng* r, unsigned mask) {
  int total = 0;
  for (int o = 0; o < OP_COUNT; ++o)
    if (mask & (1u << o)) total += fuzzOpWeight[o];
  int draw = static_cast<int>(fuzzInt(r, static_cast<size_t>(total)));
  for (int o = 0; o < OP_COUNT; ++o) {
    if (!(mask & (1u << o))) continue;
    draw -= fuzzOpWeight[o];
    if (draw < 0) return o;
  }
  return OP_RUN;
}

struct ColSpec {
  ColumnType type = ColumnType::ordinal;
  int categories = 0;      // categorical only
  double missingRate = 0.0;
};

struct ConfigSpec {
  const char* name;
  ResponseFamily family;
  std::vector<ColSpec> cols;
  size_t numChains;
  unsigned opMask;
  // > 0 builds a heteroscedastic sampler (constant-leaf gaussian only); the
  // mask must then exclude every predictor-side op, which the bridge refuses
  // under a variance forest until the routing repairs land.
  size_t numVarianceTrees = 0;
};

// One column of candidate values. Ordinal flavors: 0 identity, 1 jitter,
// 3 degenerate constant, else fresh uniform. Categorical: valid codes with
// every category present, flavor 1 injecting an out-of-range code (validating
// ops reject it), flavor 2 injecting missings.
static void fuzzFillColumn(ext_rng* r, const ColSpec& col, const double* current,
                           size_t n, int flavor, double* out) {
  const double na = std::nan("");
  if (col.type == ColumnType::categorical) {
    for (size_t i = 0; i < n; ++i)
      out[i] = static_cast<double>(fuzzInt(r, static_cast<size_t>(col.categories)));
    if (col.missingRate > 0.0 || flavor == 2)
      for (size_t i = 0; i < n; ++i)
        if (fuzzUnif(r) < (col.missingRate > 0.0 ? col.missingRate : 0.1))
          out[i] = na;
    for (int c = 0; c < col.categories && static_cast<size_t>(c) < n; ++c)
      out[c] = static_cast<double>(c);
    if (flavor == 1)
      out[fuzzInt(r, n)] = static_cast<double>(col.categories + 2);
    return;
  }
  if (flavor == 0 && current != nullptr) {
    for (size_t i = 0; i < n; ++i) out[i] = current[i];
  } else if (flavor == 1 && current != nullptr) {
    for (size_t i = 0; i < n; ++i)
      out[i] = current[i] + 0.01 * (fuzzUnif(r) - 0.5);
  } else if (flavor == 3) {
    double v = fuzzUnif(r);
    for (size_t i = 0; i < n; ++i) out[i] = v;
  } else {
    for (size_t i = 0; i < n; ++i) out[i] = fuzzUnif(r);
  }
  if (col.missingRate > 0.0)
    for (size_t i = 0; i < n; ++i)
      if (fuzzUnif(r) < col.missingRate) out[i] = na;
}

static void fuzzFillResponse(ext_rng* r, ResponseFamily fam, const double* x,
                             size_t n, size_t p, double* y) {
  for (size_t i = 0; i < n; ++i) {
    double eta = 0.0;
    for (size_t j = 0; j < p; ++j) {
      double v = x != nullptr ? x[i + j * n] : 0.0;
      if (std::isnan(v)) eta += 1.0;
      else eta += (j == 0 ? 2.0 : 1.0) * (v - 0.5);
    }
    if (fam == ResponseFamily::gaussian)
      y[i] = eta + 0.2 * (fuzzUnif(r) - 0.5);
    else
      y[i] = fuzzUnif(r) < 1.0 / (1.0 + std::exp(-eta)) ? 1.0 : 0.0;
  }
}

// The op loop, generic over the leaf model so the linear-leaf configuration
// reuses it. Returns false (and prints the seed + a trailing op trace) on the
// first invariant break.
// initialX is the dense predictor matrix the sampler was built from, or null
// for a CSC-backed sampler. The engine keeps no matrix, so the driver tracks
// the current predictors itself (as the R layer does) to seed candidate columns
// and feed the re-quantize ops (setCutPoints, setState) their raw values.
template <typename S>
static bool fuzzDrive(S& s, const ConfigSpec& spec, FuzzArena& arena,
                      ext_rng* opRng, std::uint32_t seed, int numOps,
                      const double* initialX) {
  std::vector<std::string> trace;
  char line[128];
  bool ok = true;
  int op = 0;
  auto record = [&](const char* text) { trace.push_back(text); };
  auto fail = [&](const char* what) {
    ++failures;
    printf("FAIL: fuzz [%s] seed %u op %d: %s\n", spec.name, seed, op, what);
    size_t start = trace.size() > 25 ? trace.size() - 25 : 0;
    for (size_t k = start; k < trace.size(); ++k)
      printf("  %s\n", trace[k].c_str());
    ok = false;
  };
  size_t p = s.numPredictors();
  std::vector<size_t> ordinalCols;
  for (size_t j = 0; j < p; ++j)
    if (spec.cols[j].type == ColumnType::ordinal) ordinalCols.push_back(j);

  // the driver's own copy of the current predictors, maintained in lockstep
  // with accepted mutations (empty for a CSC-backed sampler)
  std::vector<double> currentX;
  if (initialX != nullptr)
    currentX.assign(initialX, initialX + s.numObservations() * p);
  auto curCol = [&](size_t j, size_t n) -> const double* {
    return currentX.empty() ? nullptr : currentX.data() + j * n;
  };
  auto curAll = [&]() -> const double* {
    return currentX.empty() ? nullptr : currentX.data();
  };

  for (op = 0; op < numOps && ok; ++op) {
    size_t n = s.numObservations();
    int kind = fuzzPickOp(opRng, spec.opMask);
    switch (kind) {
      case OP_SET_PREDICTOR: {
        std::vector<double> cand(n * p);
        for (size_t j = 0; j < p; ++j)
          fuzzFillColumn(opRng, spec.cols[j], curCol(j, n), n,
                         static_cast<int>(fuzzInt(opRng, 5)),
                         cand.data() + j * n);
        bool force = fuzzUnif(opRng) < 0.25;
        bool upc = fuzzUnif(opRng) < 0.4;
        double* buf = arena.keep(std::move(cand));
        FuzzSnapshot<S> before = fuzzCapture(s);
        PredictorUpdateResult res = s.setPredictor(buf, force, upc);
        snprintf(line, sizeof line, "op%d setPredictor force=%d upc=%d -> %d",
                 op, force, upc, static_cast<int>(res));
        record(line);
        if (res == PredictorUpdateResult::accepted && !currentX.empty())
          currentX.assign(buf, buf + n * p);
        if (res != PredictorUpdateResult::accepted &&
            !fuzzSnapshotsEqual(before, fuzzCapture(s)))
          fail("rejected setPredictor mutated state");
        break;
      }
      case OP_UPDATE_COLUMNS: {
        size_t nc = 1 + fuzzInt(opRng, p);
        std::vector<size_t> perm(p);
        for (size_t j = 0; j < p; ++j) perm[j] = j;
        for (size_t k = 0; k < nc; ++k)
          std::swap(perm[k], perm[k + fuzzInt(opRng, p - k)]);
        std::vector<double> cols(n * nc);
        for (size_t k = 0; k < nc; ++k)
          fuzzFillColumn(opRng, spec.cols[perm[k]], curCol(perm[k], n), n,
                         static_cast<int>(fuzzInt(opRng, 5)), cols.data() + k * n);
        bool force = fuzzUnif(opRng) < 0.25;
        bool upc = fuzzUnif(opRng) < 0.4;
        FuzzSnapshot<S> before = fuzzCapture(s);
        PredictorUpdateResult res =
          s.updatePredictor(cols.data(), perm.data(), nc, force, upc);
        snprintf(line, sizeof line, "op%d updateColumns nc=%zu force=%d -> %d",
                 op, nc, force, static_cast<int>(res));
        record(line);
        if (res == PredictorUpdateResult::accepted && !currentX.empty())
          for (size_t k = 0; k < nc; ++k)
            std::memcpy(currentX.data() + perm[k] * n, cols.data() + k * n,
                        n * sizeof(double));
        if (res != PredictorUpdateResult::accepted &&
            !fuzzSnapshotsEqual(before, fuzzCapture(s)))
          fail("rejected updateColumns mutated state");
        break;
      }
      case OP_PER_OBS:
      case OP_SESSION_ABANDON: {
        size_t col = fuzzInt(opRng, p);
        int flavor = spec.cols[col].type == ColumnType::categorical
                       ? (fuzzUnif(opRng) < 0.5 ? 0 : 2)
                       : static_cast<int>(fuzzInt(opRng, 5));
        std::vector<double> nv(n);
        fuzzFillColumn(opRng, spec.cols[col], curCol(col, n), n, flavor,
                       nv.data());
        if (kind == OP_SESSION_ABANDON) {
          FuzzSnapshot<S> before = fuzzCapture(s);
          {
            std::unique_ptr<PredictorUpdateSession> session =
              s.beginPredictorUpdate(nv.data(), col);
            for (size_t i = 0; i < n; ++i)
              session->observationWouldRemainValid(i);
          }  // dropped without commit or finalize
          snprintf(line, sizeof line, "op%d sessionAbandon col=%zu", op, col);
          record(line);
          if (!fuzzSnapshotsEqual(before, fuzzCapture(s)))
            fail("abandoned session mutated state");
        } else {
          std::unique_ptr<bool[]> installed(new bool[n]);
          bool fin =
            s.updatePredictorPerObservation(nv.data(), col, installed.get());
          snprintf(line, sizeof line, "op%d perObs col=%zu fin=%d", op, col,
                   fin);
          record(line);
          if (!currentX.empty())
            for (size_t i = 0; i < n; ++i)
              if (installed[i]) currentX[col * n + i] = nv[i];
          if (!fin) fail("per-observation finalize invalid");
        }
        break;
      }
      case OP_SET_DATA: {
        size_t n2 = 24 + fuzzInt(opRng, 2 * n);
        std::vector<double> nx(n2 * p), ny(n2);
        for (size_t j = 0; j < p; ++j)
          fuzzFillColumn(opRng, spec.cols[j], nullptr, n2,
                         spec.cols[j].type == ColumnType::categorical ? 0 : 2,
                         nx.data() + j * n2);
        fuzzFillResponse(opRng, spec.family, nx.data(), n2, p, ny.data());
        double* xb = arena.keep(std::move(nx));
        double* yb = arena.keep(std::move(ny));
        s.setData(xb, yb, n2, nullptr, nullptr, nullptr, 0);
        if (!currentX.empty()) currentX.assign(xb, xb + n2 * p);
        snprintf(line, sizeof line, "op%d setData n=%zu", op, n2);
        record(line);
        break;
      }
      case OP_SET_CUTS: {
        if (ordinalCols.empty()) { record("op setCuts (skipped)"); break; }
        size_t col = ordinalCols[fuzzInt(opRng, ordinalCols.size())];
        std::uint32_t m = static_cast<std::uint32_t>(2 + fuzzInt(opRng, 6));
        std::vector<double> cuts(m);
        double v = 0.03 + 0.02 * fuzzUnif(opRng);
        for (std::uint32_t k = 0; k < m; ++k) {
          cuts[k] = v;
          v += 0.02 + 0.12 * fuzzUnif(opRng);
        }
        const double* cutPtr = cuts.data();
        s.setCutPoints(&cutPtr, &m, &col, 1, curAll());
        snprintf(line, sizeof line, "op%d setCuts col=%zu m=%u", op, col, m);
        record(line);
        break;
      }
      case OP_SET_SIGMA: {
        double sig = 0.1 + 2.0 * fuzzUnif(opRng);
        s.setSigma(sig);
        snprintf(line, sizeof line, "op%d %s", op, fuzzOpName[kind]);
        record(line);
        break;
      }
      case OP_SET_RESPONSE: {
        std::vector<double> ny(n);
        fuzzFillResponse(opRng, spec.family, curAll(), n, p, ny.data());
        double* yb = arena.keep(std::move(ny));
        s.setResponse(yb, fuzzUnif(opRng) < 0.5);
        snprintf(line, sizeof line, "op%d %s", op, fuzzOpName[kind]);
        record(line);
        break;
      }
      case OP_SET_WEIGHTS: {
        std::vector<double> nw(n);
        for (size_t i = 0; i < n; ++i) nw[i] = 0.5 + fuzzUnif(opRng);
        double* wb = arena.keep(std::move(nw));
        s.setWeights(wb);
        snprintf(line, sizeof line, "op%d %s", op, fuzzOpName[kind]);
        record(line);
        break;
      }
      case OP_SET_OFFSET: {
        std::vector<double> no(n);
        for (size_t i = 0; i < n; ++i) no[i] = 0.5 * (fuzzUnif(opRng) - 0.5);
        double* ob = arena.keep(std::move(no));
        s.setOffset(ob, fuzzUnif(opRng) < 0.5);
        snprintf(line, sizeof line, "op%d %s", op, fuzzOpName[kind]);
        record(line);
        break;
      }
      case OP_SET_TEST: {
        size_t nt = 10 + fuzzInt(opRng, 30);
        std::vector<double> xt(nt * p);
        for (size_t j = 0; j < p; ++j)
          fuzzFillColumn(opRng, spec.cols[j], nullptr, nt,
                         spec.cols[j].type == ColumnType::categorical ? 0 : 2,
                         xt.data() + j * nt);
        double* tb = arena.keep(std::move(xt));
        s.setTestPredictors(tb, nt);
        snprintf(line, sizeof line, "op%d setTest nt=%zu", op, nt);
        record(line);
        break;
      }
      case OP_RUN: {
        size_t nc = s.numChains();
        // a multi-location combiner (the multinomial softmax) reports L
        // channels per observation per sample; L is 1 everywhere else
        size_t numLocations = s.numReportedLocations();
        std::vector<double> sig(nc), tf(nc * n * numLocations);
        Results r;
        r.numReportedLocations = numLocations;
        r.sigma = sig.data();
        r.trainingFits = tf.data();
        s.run(0, 1, r);
        snprintf(line, sizeof line, "op%d %s", op, fuzzOpName[kind]);
        record(line);
        bool finite = true;
        for (double v : tf) finite &= std::isfinite(v);
        for (double v : sig) finite &= std::isfinite(v);
        if (!finite) fail("run produced non-finite draws");
        if (spec.family == ResponseFamily::gaussian)
          for (double v : sig)
            if (!(v > 0.0)) fail("run produced non-positive sigma");
        break;
      }
      case OP_STATE: {
        // Own-state round trip through serialization: a state a live sampler
        // reaches must restore into itself. All the mutation edges that once
        // broke this are fixed at their sources.
        SamplerStateData st;
        s.getState(st);
        record("op state");
        if (!s.setState(st, curAll())) {
          fail("setState refused own state");
        } else {
          SamplerStateData st2;
          s.getState(st2);
          if (!statesAgree(st, st2)) fail("state round trip disagrees");
        }
        break;
      }
      default: break;
    }
    if (ok) {
      const char* v = fuzzInvariantViolation(s);
      if (v) fail(v);
    }
  }
  return ok;
}

static const double fuzzRawScale = 0.37804942330213542;  // qchisq(0.1, 3) / 3

static void fuzzRunConstant(const ConfigSpec& spec, std::uint32_t seed,
                            int numOps) {
  FuzzArena arena;
  size_t n0 = 160, p = spec.cols.size();
  ext_rng* opRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(opRng, seed * 2u + 7u);

  std::vector<ColumnType> types(p);
  std::vector<double> x(n0 * p), y(n0);
  for (size_t j = 0; j < p; ++j) {
    types[j] = spec.cols[j].type;
    fuzzFillColumn(opRng, spec.cols[j], nullptr, n0,
                   spec.cols[j].type == ColumnType::categorical ? 0 : 2,
                   x.data() + j * n0);
  }
  fuzzFillResponse(opRng, spec.family, x.data(), n0, p, y.data());
  double* xb = arena.keep(std::move(x));
  double* yb = arena.keep(std::move(y));

  SamplerOptions options;
  options.numTrees = 15;
  options.numChains = spec.numChains;
  options.numVarianceTrees = spec.numVarianceTrees;
  options.predictors.columnTypes = types.data();
  if (spec.family != ResponseFamily::gaussian) options.nodeScale = 3.0;
  std::vector<ext_rng*> rngs(spec.numChains);
  for (size_t c = 0; c < spec.numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], seed * 2u + 1u + static_cast<std::uint32_t>(c));
  }
  double rawScale =
    spec.family == ResponseFamily::gaussian ? fuzzRawScale : 1.0;
  Sampler<ConstantGaussianLeaf> s(xb, yb, n0, p, nullptr, nullptr, spec.family,
                                  1.0, 3.0, rawScale, options, rngs.data());
  Results empty;
  s.run(12, 0, empty);
  fuzzDrive(s, spec, arena, opRng, seed, numOps, xb);

  for (size_t c = 0; c < spec.numChains; ++c) ext_rng_destroy(rngs[c]);
  ext_rng_destroy(opRng);
}

// The multi-forest op surface: the transactional predictor paths, the forced
// refresh, setCutPoints, run and the state round trip. setData, setResponse,
// setWeights and setOffset are refused for a multi-forest sampler at the
// bridge (test-multi-forest-seam.R), and the per-observation session joins at
// S2 of docs/plans/multiforest-predictor-mutation.md, when its cell guard
// widens past forest 0.
static const unsigned fuzzMultiForestMask =
  (1u << OP_SET_PREDICTOR) | (1u << OP_UPDATE_COLUMNS) | (1u << OP_SET_CUTS) |
  (1u << OP_RUN) | (1u << OP_STATE);

// A two-forest BCF sampler (docs/design/bcf.md) over the same op loop: the
// shape whose forest 1 the widened revalidation must re-route, and which the
// snapshot above must see. tau carries a moderator subset, so its trees split
// on a strict subset of the columns and the j-split pruning of the rebuild is
// exercised rather than vacuous.
static void fuzzRunBCF(const ConfigSpec& spec, std::uint32_t seed, int numOps) {
  FuzzArena arena;
  size_t n0 = 160, p = spec.cols.size();
  ext_rng* opRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(opRng, seed * 2u + 7u);

  std::vector<ColumnType> types(p);
  std::vector<double> x(n0 * p), y(n0), z(n0);
  for (size_t j = 0; j < p; ++j) {
    types[j] = spec.cols[j].type;
    fuzzFillColumn(opRng, spec.cols[j], nullptr, n0,
                   spec.cols[j].type == ColumnType::categorical ? 0 : 2,
                   x.data() + j * n0);
  }
  fuzzFillResponse(opRng, spec.family, x.data(), n0, p, y.data());
  for (size_t i = 0; i < n0; ++i) {
    z[i] = fuzzUnif(opRng) < 0.5 ? 0.0 : 1.0;
    y[i] += z[i] * (1.0 + 2.0 * x[i]);
  }
  double* xb = arena.keep(std::move(x));
  double* yb = arena.keep(std::move(y));
  double* zb = arena.keep(std::move(z));

  SamplerOptions options;
  options.numChains = spec.numChains;
  options.predictors.columnTypes = types.data();
  std::vector<size_t> moderators = {0, 1};
  BCFSpec bcf;
  bcf.mu.numTrees = 15;
  bcf.tau.numTrees = 10;
  bcf.tau.base = 0.25;
  bcf.tau.power = 3.0;
  bcf.tau.columns = moderators.data();
  bcf.tau.numColumns = moderators.size();
  bcf.z = zb;

  std::vector<ext_rng*> rngs(spec.numChains);
  for (size_t c = 0; c < spec.numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], seed * 2u + 1u + static_cast<std::uint32_t>(c));
  }
  Sampler<ConstantGaussianLeaf> s(xb, yb, n0, p, nullptr, nullptr, 1.0, 3.0,
                                  fuzzRawScale, options, bcf, rngs.data());
  Results empty;
  s.run(12, 0, empty);
  fuzzDrive(s, spec, arena, opRng, seed, numOps, xb);

  for (size_t c = 0; c < spec.numChains; ++c) ext_rng_destroy(rngs[c]);
  ext_rng_destroy(opRng);
}

// A K-forest multinomial (softmax) sampler (docs/design/multinomial.md): the
// other shape the widened revalidation lifts, and the one with the largest
// j-splitting tree count. Its response is the borrowed one-hot count matrix,
// so the response-side ops stay out for that reason as well.
static void fuzzRunMultinomial(const ConfigSpec& spec, std::uint32_t seed,
                               int numOps) {
  FuzzArena arena;
  size_t n0 = 160, p = spec.cols.size(), K = 3;
  ext_rng* opRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(opRng, seed * 2u + 7u);

  std::vector<ColumnType> types(p);
  std::vector<double> x(n0 * p);
  for (size_t j = 0; j < p; ++j) {
    types[j] = spec.cols[j].type;
    fuzzFillColumn(opRng, spec.cols[j], nullptr, n0,
                   spec.cols[j].type == ColumnType::categorical ? 0 : 2,
                   x.data() + j * n0);
  }
  // a covariate-dependent label, one-hot into the category-major count matrix
  std::vector<int> counts(n0 * K, 0), trials(n0, 1);
  for (size_t i = 0; i < n0; ++i) {
    double v = x[i];
    if (std::isnan(v)) v = 0.5;
    size_t k = fuzzUnif(opRng) < 0.5 + 0.4 * (v - 0.5)
                 ? 0
                 : (fuzzUnif(opRng) < 0.5 ? 1 : 2);
    counts[k * n0 + i] = 1;
  }
  double* xb = arena.keep(std::move(x));

  SamplerOptions options;
  options.numChains = spec.numChains;
  options.predictors.columnTypes = types.data();
  MultinomialSpec mn;
  mn.numCategories = K;
  mn.counts = counts.data();
  mn.trials = trials.data();
  mn.forest.numTrees = 12;

  std::vector<ext_rng*> rngs(spec.numChains);
  for (size_t c = 0; c < spec.numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], seed * 2u + 1u + static_cast<std::uint32_t>(c));
  }
  Sampler<ConstantGaussianLeaf> s(xb, n0, p, options, mn, rngs.data());
  Results empty;
  s.run(12, 0, empty);
  fuzzDrive(s, spec, arena, opRng, seed, numOps, xb);

  for (size_t c = 0; c < spec.numChains; ++c) ext_rng_destroy(rngs[c]);
  ext_rng_destroy(opRng);
}

static void fuzzRunLinear(std::uint32_t seed, int numOps) {
  FuzzArena arena;
  size_t n0 = 160, p = 3;
  ConfigSpec spec{"linear", ResponseFamily::gaussian,
                  {ColSpec{}, ColSpec{}, ColSpec{}}, 1,
                  (1u << OP_COUNT) - 1};
  ext_rng* opRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(opRng, seed * 2u + 7u);

  std::vector<double> x(n0 * p), y(n0);
  for (size_t j = 0; j < p; ++j)
    fuzzFillColumn(opRng, spec.cols[j], nullptr, n0, 2, x.data() + j * n0);
  fuzzFillResponse(opRng, spec.family, x.data(), n0, p, y.data());
  double* xb = arena.keep(std::move(x));
  double* yb = arena.keep(std::move(y));

  std::vector<size_t> covariates = {2};
  SamplerOptions options;
  options.numTrees = 15;
  options.leafCovariateColumns = covariates.data();
  options.numLeafCovariates = 1;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, seed * 2u + 1u);
  Sampler<LinearGaussianLeaf> s(xb, yb, n0, p, nullptr, nullptr,
                                ResponseFamily::gaussian, 1.0, 3.0, fuzzRawScale,
                                options, &rng);
  Results empty;
  s.run(12, 0, empty);
  fuzzDrive(s, spec, arena, opRng, seed, numOps, xb);

  ext_rng_destroy(rng);
  ext_rng_destroy(opRng);
}

static void fuzzRunSparse(std::uint32_t seed, int numOps) {
  FuzzArena arena;
  const size_t n0 = 200;
  ConfigSpec spec{"sparse", ResponseFamily::gaussian,
                  {ColSpec{}, ColSpec{}, ColSpec{}, ColSpec{}}, 1,
                  (1u << OP_SET_SIGMA) | (1u << OP_SET_RESPONSE) |
                    (1u << OP_SET_WEIGHTS) | (1u << OP_RUN) | (1u << OP_STATE)};
  ext_rng* opRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(opRng, seed * 2u + 7u);

  rngState = 0x9E3779B97F4A7C15ull ^ static_cast<uint64_t>(seed);  // CscFixture
  CscFixture fixture;
  fixture.build(n0, {0.15, 0.15, 0.15, 0.15}, 4);
  std::vector<double> y(n0);
  fuzzFillResponse(opRng, spec.family, fixture.dense.data(), n0, fixture.p,
                   y.data());
  double* yb = arena.keep(std::move(y));

  SamplerOptions options;
  options.numTrees = 15;
  options.predictors.cscColumnPointers = fixture.pointers.data();
  options.predictors.cscRowIndices = fixture.rows.data();
  options.predictors.cscValues = fixture.values.data();
  options.predictors.columnSources = fixture.allCscSources.data();
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, seed * 2u + 1u);
  Sampler<ConstantGaussianLeaf> s(nullptr, yb, n0, fixture.p, nullptr, nullptr,
                                  ResponseFamily::gaussian, 1.0, 3.0,
                                  fuzzRawScale, options, &rng);
  Results empty;
  s.run(12, 0, empty);
  // CSC-backed: the engine re-quantizes from its retained slices, so the driver
  // tracks no dense matrix
  fuzzDrive(s, spec, arena, opRng, seed, numOps, nullptr);

  ext_rng_destroy(rng);
  ext_rng_destroy(opRng);
}

// F4, the structural gate on the snapshot (docs/plans/multiforest-predictor-
// mutation.md): table-driven, one row per captured family, each row perturbing
// exactly that family and requiring the comparison to go false. A family with
// no row is a family the snapshot does not cover - which is the failure mode
// the rebuild exists to close, and which no amount of green fuzzing detects.
using FuzzSnap = FuzzSnapshot<Sampler<ConstantGaussianLeaf>>;

struct FuzzFamily {
  const char* name;
  std::function<void(FuzzSnap&)> perturb;
};

// pre-order index of the first leaf of a flattened tree; every tree has one
static size_t fuzzFirstFlatLeaf(const std::vector<FlatNode>& tree) {
  for (size_t i = 0; i < tree.size(); ++i)
    if (tree[i].variable == invalidVariable) return i;
  return 0;
}

// arena index of the first internal node of a live tree, or 0 for a stump
static size_t fuzzFirstInternal(const std::vector<FuzzNode>& nodes) {
  for (size_t i = 0; i < nodes.size(); ++i)
    if (nodes[i].leftChild != invalidNode) return i;
  return 0;
}

static void fuzzCheckFamilies(FuzzSnap& base,
                              const std::vector<FuzzFamily>& families) {
  char line[160];
  for (const FuzzFamily& family : families) {
    FuzzSnap moved = base;
    family.perturb(moved);
    snprintf(line, sizeof line,
             "the fuzz snapshot sees a change in %s", family.name);
    check(!fuzzSnapshotsEqual(base, moved), line);
  }
}

static void testSnapshotCoversEveryFamily() {
  const size_t n = 120, p = 3;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 0.0 : 1.0;
    y[i] = 4.0 * (x[i] - 0.5) + 2.0 * x[i + n] +
           z[i] * (1.0 + 2.0 * x[i + 2 * n]) + 0.2 * (runif01() - 0.5);
  }

  // a two-chain BCF: the mean/tau, per-chain and persisted families
  {
    SamplerOptions options;
    options.numChains = 2;
    BCFSpec spec;
    spec.mu.numTrees = 10;
    spec.tau.numTrees = 8;
    spec.z = z.data();
    ext_rng* rngs[2];
    for (size_t c = 0; c < 2; ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], 5150u + static_cast<std::uint32_t>(c));
    }
    Sampler<ConstantGaussianLeaf> s(x.data(), y.data(), n, p, nullptr, nullptr,
                                    1.0, 3.0, fuzzRawScale, options, spec, rngs);
    Results empty;
    s.run(60, 0, empty);
    FuzzSnap base = fuzzCapture(s);
    check(base.geom.trees.size() == 2 && base.geom.trees[0].size() == 2,
          "the snapshot captured both chains and both BCF forests");
    check(fuzzFirstInternal(base.geom.trees[0][1][0].nodes) != 0 ||
            base.geom.trees[0][1][0].nodes.size() > 1,
          "the BCF snapshot's tau trees carry a split to perturb");

    fuzzCheckFamilies(base, {
      {"a store code", [](FuzzSnap& t) { t.codes[0] ^= 1; }},
      {"a cut point", [](FuzzSnap& t) { t.state.cutPoints[0][0] += 1.0; }},
      {"the saved-tree write position",
       [](FuzzSnap& t) { ++t.state.currentSampleNum; }},
      {"sigma", [](FuzzSnap& t) { t.state.chains[0].sigma += 1.0; }},
      {"the BCF glue", [](FuzzSnap& t) { t.state.chains[0].a += 1.0; }},
      {"a persisted mean leaf", [](FuzzSnap& t) {
         std::vector<FlatNode>& tree(t.state.chains[0].forests[0].trees[0]);
         tree[fuzzFirstFlatLeaf(tree)].value += 1.0;
       }},
      {"a persisted tau leaf", [](FuzzSnap& t) {
         std::vector<FlatNode>& tree(t.state.chains[0].forests[1].trees[0]);
         tree[fuzzFirstFlatLeaf(tree)].value += 1.0;
       }},
      {"a mean split rule", [](FuzzSnap& t) {
         std::vector<FuzzNode>& nodes(t.geom.trees[0][0][0].nodes);
         ++nodes[fuzzFirstInternal(nodes)].variable;
       }},
      {"a tau split rule", [](FuzzSnap& t) {
         std::vector<FuzzNode>& nodes(t.geom.trees[0][1][0].nodes);
         nodes[fuzzFirstInternal(nodes)].bits ^= 1;
       }},
      {"a mean partition range",
       [](FuzzSnap& t) { ++t.geom.trees[0][0][0].nodes[0].end; }},
      {"a tau partition range",
       [](FuzzSnap& t) { ++t.geom.trees[0][1][0].nodes[0].end; }},
      {"a mean leaf assignment",
       [](FuzzSnap& t) { ++t.geom.trees[0][0][0].leafOwner[0]; }},
      {"a tau leaf assignment",
       [](FuzzSnap& t) { ++t.geom.trees[0][1][0].leafOwner[0]; }},
      {"a pooled mask word",
       [](FuzzSnap& t) { t.geom.trees[0][0][0].maskPool.push_back(1ull); }},
      {"mean tree fits", [](FuzzSnap& t) { t.geom.fits[0][0][0] += 1.0; }},
      {"tau tree fits", [](FuzzSnap& t) { t.geom.fits[0][1][0] += 1.0; }},
      {"mean totalFits", [](FuzzSnap& t) { t.geom.totals[0][0][0] += 1.0; }},
      {"tau totalFits", [](FuzzSnap& t) { t.geom.totals[0][1][0] += 1.0; }},
      {"chain 1's tau fits", [](FuzzSnap& t) { t.geom.fits[1][1][0] += 1.0; }},
      {"chain 1's tau partition",
       [](FuzzSnap& t) { ++t.geom.trees[1][1][0].leafOwner[0]; }},
    });
    for (size_t c = 0; c < 2; ++c) ext_rng_destroy(rngs[c]);
  }

  // a heteroscedastic sampler: the variance-forest families, which live
  // outside forests_ and which the old snapshot could not see at all
  {
    SamplerOptions options;
    options.numTrees = 10;
    options.numVarianceTrees = 8;
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 5160u);
    Sampler<ConstantGaussianLeaf> s(x.data(), y.data(), n, p, nullptr, nullptr,
                                    ResponseFamily::gaussian, 1.0, 3.0,
                                    fuzzRawScale, options, &rng);
    Results empty;
    s.run(60, 0, empty);
    FuzzSnap base = fuzzCapture(s);
    check(base.geom.varianceTrees[0].size() == 8,
          "the snapshot captured every variance tree");

    fuzzCheckFamilies(base, {
      {"a persisted variance leaf", [](FuzzSnap& t) {
         std::vector<FlatNode>& tree(t.state.chains[0].varianceTrees[0]);
         tree[fuzzFirstFlatLeaf(tree)].value += 1.0;
       }},
      {"a variance split rule", [](FuzzSnap& t) {
         std::vector<FuzzNode>& nodes(t.geom.varianceTrees[0][0].nodes);
         nodes[fuzzFirstInternal(nodes)].bits ^= 1;
       }},
      {"a variance partition range",
       [](FuzzSnap& t) { ++t.geom.varianceTrees[0][0].nodes[0].end; }},
      {"a variance leaf assignment",
       [](FuzzSnap& t) { ++t.geom.varianceTrees[0][0].leafOwner[0]; }},
      {"a variance factor",
       [](FuzzSnap& t) { t.geom.varianceFactors[0][0] += 1.0; }},
      {"the combined variance",
       [](FuzzSnap& t) { t.geom.varianceFits[0][0] += 1.0; }},
    });
    ext_rng_destroy(rng);
  }

  printf("ok: fuzz snapshot family coverage\n");
}

static void testMutationFuzzer(int numSeeds) {
  ColSpec ord;
  ColSpec cat4{ColumnType::categorical, 4, 0.0};
  ColSpec ordMiss{ColumnType::ordinal, 0, 0.25};
  ColSpec cat3Miss{ColumnType::categorical, 3, 0.2};
  std::vector<ConfigSpec> configs = {
    {"gaussian", ResponseFamily::gaussian, {ord, ord, ord}, 1,
     (1u << OP_COUNT) - 1},
    {"categorical", ResponseFamily::gaussian, {cat4, ord, ord}, 1,
     (1u << OP_COUNT) - 1},
    {"missing", ResponseFamily::gaussian, {ordMiss, cat3Miss, ord}, 1,
     (1u << OP_COUNT) - 1},
    {"multichain", ResponseFamily::gaussian, {ord, ord, ord}, 2,
     (1u << OP_COUNT) - 1},
    {"probit", ResponseFamily::probit, {ord, ord, ord}, 1,
     (1u << OP_SET_PREDICTOR) | (1u << OP_UPDATE_COLUMNS) | (1u << OP_PER_OBS) |
       (1u << OP_SESSION_ABANDON) | (1u << OP_SET_DATA) | (1u << OP_SET_CUTS) |
       (1u << OP_SET_TEST) | (1u << OP_RUN) | (1u << OP_STATE)},
    {"logistic", ResponseFamily::logistic, {ord, ord, ord}, 1,
     (1u << OP_SET_PREDICTOR) | (1u << OP_UPDATE_COLUMNS) | (1u << OP_PER_OBS) |
       (1u << OP_SESSION_ABANDON) | (1u << OP_SET_DATA) | (1u << OP_SET_CUTS) |
       (1u << OP_SET_TEST) | (1u << OP_RUN) | (1u << OP_STATE)},
    // setCuts joined the safe surface once setCutPoints re-routed the variance
    // forest through forceRefreshTrees, and setData once applyNewData resized
    // its seven n-sized allocations and re-routed it - the op the fuzzer runs
    // at a fresh count every time. The transactional predictor ops and
    // setSigma stay out for the whole arc: they are refused at the bridge
    // (setSigma is additionally an engine no-op under a variance forest).
    {"heteroscedastic", ResponseFamily::gaussian, {ord, ord, ord}, 1,
     (1u << OP_SET_RESPONSE) | (1u << OP_SET_WEIGHTS) | (1u << OP_SET_OFFSET) |
       (1u << OP_SET_TEST) | (1u << OP_RUN) | (1u << OP_STATE) |
       (1u << OP_SET_CUTS) | (1u << OP_SET_DATA),
     8},
  };
  // the multi-forest shapes, on their own runners: one continuous and one
  // categorical design each, so the transactional paths meet both column types
  ConfigSpec bcf{"bcf", ResponseFamily::gaussian, {ord, ord, ord}, 1,
                 fuzzMultiForestMask};
  ConfigSpec bcfCat{"bcf-categorical", ResponseFamily::gaussian,
                    {ord, cat4, ord}, 2, fuzzMultiForestMask};
  ConfigSpec multinomial{"multinomial", ResponseFamily::logistic,
                         {ord, ord, ord}, 1, fuzzMultiForestMask};
  ConfigSpec multinomialMiss{"multinomial-missing", ResponseFamily::logistic,
                             {ordMiss, cat3Miss, ord}, 1, fuzzMultiForestMask};

  const int numOps = 40;
  for (int sd = 0; sd < numSeeds; ++sd) {
    std::uint32_t seed = 1u + static_cast<std::uint32_t>(sd);
    for (const ConfigSpec& spec : configs)
      fuzzRunConstant(spec, seed, numOps);
    fuzzRunBCF(bcf, seed, numOps);
    fuzzRunBCF(bcfCat, seed, numOps);
    fuzzRunMultinomial(multinomial, seed, numOps);
    fuzzRunMultinomial(multinomialMiss, seed, numOps);
    fuzzRunLinear(seed, numOps);
    fuzzRunSparse(seed, numOps);
  }
  printf("ok: mutation fuzzer (%d seeds)\n", numSeeds);
}

void runFuzzTests(int numSeeds) {
  testSnapshotCoversEveryFamily();
  testMutationFuzzer(numSeeds);
}
