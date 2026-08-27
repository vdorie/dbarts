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
  // the count matrix and its trials are borrowed int buffers, on the same
  // lifetime rule as the double ones
  std::vector<std::unique_ptr<std::vector<int>>> intBufs;
  double* keep(std::vector<double>&& v) {
    bufs.push_back(std::make_unique<std::vector<double>>(std::move(v)));
    return bufs.back()->data();
  }
  int* keep(std::vector<int>&& v) {
    intBufs.push_back(std::make_unique<std::vector<int>>(std::move(v)));
    return intBufs.back()->data();
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

// The state a rejected or abandoned mutation must leave untouched. An earlier
// fingerprint carried codes, cuts, sigma and Chain::treeFits(), which is
// forest 0 alone: a rollback that restored mu and left tau routed by the
// proposal passed it green. It is rebuilt here so that omission stops being
// expressible: the persisted state entire, plus every LIVE structure the
// persisted state does not carry, captured by INDEX LOOPS with no forest
// literal anywhere - a new forest class is covered on arrival.
// testSnapshotCoversEveryFamily is the gate that each captured family, when
// perturbed alone, makes the comparison false.

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
// its members. Measured on the shipped rollback, across both the old
// forest-0-only revalidation and the widened per-forest one: every other
// family below - codes, cuts, the whole persisted state, split
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
        trees[t] = fuzzCaptureTree(ch.treeInForest(f, t), n);
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
        fuzzCaptureTree(ch.varianceTree(j), n));
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
         a.state.recordedDraws == b.state.recordedDraws &&
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
static const char* fuzzInvariantViolation(S& s, const double* z) {
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
        const Tree& tree(ch.treeInForest(f, t));
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
    // I2, amplitude parts identity: the combined location the accessor reports
    // IS sum_f m_f(i) * forest f's totals, on the response scale. Selection is
    // EXACT rather than a min over the amplitude block - forest 0 carries the
    // synthesized all-ones basis, whose multiplier is its single amplitude, an
    // indicator forest the complementary 0/1 pair, b0 on control and b1 on
    // treated - so a z-indexing defect that swapped the two shows here, where
    // a min over candidates would absorb it. Accumulating from the LAST forest
    // DOWN matches combinedFits, so fma contraction lands on the same product
    // on both sides. It discriminates a combination gone stale against the
    // forests and amplitudes it summarizes - a cached blend, a dropped forest
    // - and is blind by construction to everything UPSTREAM of them. The same
    // identity is pinned at a fixed scenario in test-fits-without-offset.R;
    // new here is holding it after every op of a random mutation stream.
    if (z != nullptr && s.totalAmplitudes() > 0 &&
        s.numReportedLocations() == 1) {
      size_t numForests = ch.numForests(), last = numForests - 1;
      std::vector<double> amplitude(s.totalAmplitudes()), combined(n);
      std::vector<size_t> block(numForests + 1, 0);
      for (size_t f = 0; f < numForests; ++f)
        block[f + 1] = block[f] + s.numForestAmplitudes(f);
      if (!s.amplitudes(c, amplitude.data())) return "amplitudes refused";
      if (!s.fitsWithoutOffset(c, combined.data()))
        return "fitsWithoutOffset refused a single-location coupling";
      ForestCalibration calibration = s.forestCalibration(c, 0);
      auto multiplier = [&](size_t f, size_t i) {
        return s.numForestAmplitudes(f) == 1
                 ? amplitude[block[f]]
                 : amplitude[block[f] + (z[i] != 0.0 ? 1 : 0)];
      };
      for (size_t i = 0; i < n; ++i) {
        double location = multiplier(last, i) * ch.totalFitsInForest(last)[i];
        for (size_t f = last; f-- > 0;)
          location += multiplier(f, i) * ch.totalFitsInForest(f)[i];
        double expected =
          calibration.responseScale * location + calibration.responseShift;
        if (std::fabs(combined[i] - expected) >
            1e-9 * (1.0 + std::fabs(expected)))
          return "fitsWithoutOffset != the amplitude blend of forest totals";
      }
    }
    std::vector<int32_t> owner;
    for (size_t f = 0; f < ch.numForests(); ++f)
      for (size_t t = 0; t < ch.numTreesInForest(f); ++t)
        if (!fuzzTreeRoutesCorrectly(ch.treeInForest(f, t),
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
  OP_SET_OFFSET, OP_SET_TEST, OP_RUN, OP_STATE, OP_GROW, OP_SET_COUNTS,
  OP_SET_CATEGORY_OFFSET, OP_COUNT
};

static const char* const fuzzOpName[OP_COUNT] = {
  "setPredictor", "updateColumns", "perObs", "sessionAbandon", "setData",
  "setCuts", "setSigma", "setResponse", "setWeights", "setOffset", "setTest",
  "run", "state", "grow", "setCounts", "setCategoryOffset"};

static const int fuzzOpWeight[OP_COUNT] = {5, 5, 4, 3, 2, 3, 2, 2, 2, 2, 2,
                                           3, 2, 2, 3, 3};

// The multinomial response-side ops. fuzzPickOp divides a running total over
// the SET bits, so admitting these into a config that does not want them
// re-deals its whole op sequence and silently discards the sequences it
// explores today; single-forest configs name the complement, not all ones.
static const unsigned fuzzCountsOps =
  (1u << OP_SET_COUNTS) | (1u << OP_SET_CATEGORY_OFFSET);
static const unsigned fuzzAllSingleForestOps =
  ((1u << OP_COUNT) - 1) & ~fuzzCountsOps;

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
  // > 0 builds a heteroscedastic sampler (constant-leaf gaussian only)
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
// treatmentZ is the 0/1 indicator an amplitude coupling's forest 1 basis was
// synthesized from, or null; it selects I2's exact multiplier and nothing else.
template <typename S>
static bool fuzzDrive(S& s, const ConfigSpec& spec, FuzzArena& arena,
                      ext_rng* opRng, std::uint32_t seed, int numOps,
                      const double* initialX,
                      const double* treatmentZ = nullptr) {
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
  // the borrowed multinomial response-side buffers currently installed: the
  // counts to re-install for the self-swap flavor, the offset so I1 can
  // recompute the reported softmax against the same shift the combiner blends
  const int* liveCounts = nullptr;
  const int* liveTrials = nullptr;
  const double* liveCategoryOffset = nullptr;

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
        // I3, liveness: every reported draw is finite. It discriminates the
        // reporting slab that was never sized or never written - the coarsest
        // of the four, and the only one that fires before a value is
        // comparable to anything at all.
        if (!finite) fail("run produced non-finite draws");
        if (spec.family == ResponseFamily::gaussian)
          for (double v : sig)
            if (!(v > 0.0)) fail("run produced non-positive sigma");
        // I1, multinomial parts identity: the reported channel IS softmax over
        // the live per-category totals plus the installed category offset,
        // mirroring softmaxLocationMajor's max-subtraction so the two agree to
        // a few ulps rather than to exp()'s own absolute error. It
        // discriminates a reporting blend gone stale against the forests it
        // summarizes - a per-column refresh, a cached slab, a dropped offset -
        // none of which the simplex half alone would see; it is blind by
        // construction to everything UPSTREAM of totalFits, so a wrong draw
        // reported faithfully passes. The same identity is pinned at fixed
        // scenarios in test-multinomial-counts-mutation.R and
        // test-multinomial-category-offset.R; new here is holding it across a
        // random 40-op stream, multi-chain, under ASAN. Sampler::run slabs
        // trainingFits chain-major at c * numSamples * n * L, location-major
        // within a sample, so the chain stride is n * L at numSamples 1.
        std::vector<const double*> total(numLocations);
        for (size_t c = 0; ok && numLocations > 1 && c < nc; ++c) {
          const auto& ch(s.chain(c));
          for (size_t k = 0; k < numLocations; ++k)
            total[k] = ch.totalFitsInForest(k).data();
          const double* rep = tf.data() + c * n * numLocations;
          for (size_t i = 0; ok && i < n; ++i) {
            auto raw = [&](size_t k) {
              return total[k][i] + (liveCategoryOffset != nullptr
                                      ? liveCategoryOffset[k * n + i] : 0.0);
            };
            double maxRaw = raw(0), sumExp = 0.0, rowSum = 0.0;
            for (size_t k = 1; k < numLocations; ++k)
              maxRaw = std::max(maxRaw, raw(k));
            for (size_t k = 0; k < numLocations; ++k)
              sumExp += std::exp(raw(k) - maxRaw);
            for (size_t k = 0; ok && k < numLocations; ++k) {
              double q = std::exp(raw(k) - maxRaw) / sumExp;
              double got = rep[k * n + i];
              rowSum += got;
              if (!(got >= 0.0 && got <= 1.0))
                fail("reported channel outside [0, 1]");
              else if (std::fabs(got - q) > 1e-12 * (1.0 + std::fabs(q)))
                fail("reported channel != softmax of the category totals");
            }
            if (ok && std::fabs(rowSum - 1.0) > 1e-12)
              fail("reported channel row does not sum to 1");
          }
        }
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
      case OP_GROW: {
        // The warm-start initializer, run on live state. The generic invariant
        // below is the cheap oracle for a candidate that would empty a child -
        // it reports "empty leaf" - and the state round trip is the cheapest
        // one for the flat encoding of the rules grow places, categorical
        // masks included, on the configurations that carry them.
        size_t numSweeps = 1 + fuzzInt(opRng, 2);
        s.growFromRoot(numSweeps);
        snprintf(line, sizeof line, "op%d grow sweeps=%zu", op, numSweeps);
        record(line);
        // Before the round trip, not folded into the invariant call after
        // the switch: setState below re-derives every fit from the
        // deserialized trees, so it would silently heal a forest whose
        // derived fits are left inconsistent with the trees growFromRoot
        // just rebuilt - the defect class this line exists to catch. A
        // stale leafOf index from the same kind of mid-grow slip never
        // reaches this check: rollTreeResidual's own assertion aborts the
        // process first, and only when assertions are live (no NDEBUG).
        if (const char* v = fuzzInvariantViolation(s, treatmentZ)) fail(v);
        SamplerStateData grown;
        s.getState(grown);
        if (!s.setState(grown, curAll())) {
          fail("setState refused a grown state");
        } else {
          SamplerStateData reGrown;
          s.getState(reGrown);
          if (!statesAgree(grown, reGrown))
            fail("grown state round trip disagrees");
        }
        break;
      }
      // The count-matrix swap: the multinomial response-side channel, and the
      // only route by which Sampler::setCounts' per-chain fan-out is reached
      // from a test. Only VALID counts are generated - non-negative small
      // whole cells, every row's trials the recomputed sum and at least 1 -
      // because the contract is host-enforced and the engine is a two-pointer
      // swap; the hostile-input arm belongs at the bridge, which probes it.
      // n and K are read from the LIVE sampler at op time, so a later widening
      // of the multi-forest mask admitting setData could not meet a
      // stale-sized buffer here. NEITHER new op has an oracle for its own
      // EFFECT: the draws move and the fuzz compares draws to nothing, so
      // these are REACHABILITY ops whose oracles are I1, the ASAN leg, and
      // I4 - the counts ride no wire block, so a later OP_STATE in the same
      // stream is the state round trip across a swap.
      case OP_SET_COUNTS: {
        if (!s.supportsCountsMutation()) {
          record("op setCounts (skipped)");
          break;
        }
        size_t K = s.numReportedLocations();
        // flavors: 0 unit one-hot (the creation shape), 1 grouped multi-trial,
        // 2 a self-swap of the resident pointers, which must be BITWISE inert -
        // near-vacuous against today's plain assignment, and the only tripwire
        // that would survive a later count-derived cache
        int flavor = static_cast<int>(
          fuzzInt(opRng, liveCounts == nullptr ? 2 : 3));
        if (flavor == 2) {
          FuzzSnapshot<S> before = fuzzCapture(s);
          bool swapped = s.setCounts(liveCounts, liveTrials);
          record("op setCounts self-swap");
          if (!swapped) fail("setCounts refused on a counts-owning coupling");
          else if (!fuzzSnapshotsEqual(before, fuzzCapture(s)))
            fail("self setCounts mutated state");
          break;
        }
        std::vector<int> counts(n * K, 0), trials(n, 0);
        for (size_t i = 0; i < n; ++i) {
          if (flavor == 0) {
            counts[fuzzInt(opRng, K) * n + i] = 1;
          } else {
            for (size_t k = 0; k < K; ++k)
              counts[k * n + i] = static_cast<int>(fuzzInt(opRng, 3));
          }
          int rowSum = 0;
          for (size_t k = 0; k < K; ++k) rowSum += counts[k * n + i];
          if (rowSum == 0) {
            counts[fuzzInt(opRng, K) * n + i] = 1;
            rowSum = 1;
          }
          trials[i] = rowSum;
        }
        liveCounts = arena.keep(std::move(counts));
        liveTrials = arena.keep(std::move(trials));
        if (!s.setCounts(liveCounts, liveTrials))
          fail("setCounts refused on a counts-owning coupling");
        snprintf(line, sizeof line, "op%d setCounts flavor=%d", op, flavor);
        record(line);
        break;
      }
      // The n x K category offset, installed or cleared. The driver keeps the
      // live pointer because I1 has to recompute the reported softmax against
      // the same shift the combiner blends. OP_SET_CATEGORY_TEST_OFFSET has no
      // op here: these configs never set test predictors (OP_SET_TEST is
      // outside fuzzMultiForestMask), so numTestObservations is 0 and a test
      // offset would install a pointer nothing reads.
      case OP_SET_CATEGORY_OFFSET: {
        if (!s.supportsCountsMutation()) {
          record("op setCategoryOffset (skipped)");
          break;
        }
        size_t K = s.numReportedLocations();
        bool clear = fuzzUnif(opRng) < 0.3;
        const double* offset = nullptr;
        if (!clear) {
          std::vector<double> values(n * K);
          for (size_t m = 0; m < n * K; ++m)
            values[m] = 0.8 * (fuzzUnif(opRng) - 0.5);
          offset = arena.keep(std::move(values));
        }
        if (!s.setCategoryOffset(offset))
          fail("setCategoryOffset refused on a counts-owning coupling");
        liveCategoryOffset = offset;
        snprintf(line, sizeof line, "op%d setCategoryOffset clear=%d", op,
                 clear);
        record(line);
        break;
      }
      default: break;
    }
    if (ok) {
      const char* v = fuzzInvariantViolation(s, treatmentZ);
      if (v) fail(v);
    }
  }
  // A deterministic fingerprint of the whole op stream: every op path records
  // a line carrying its index, its drawn sizes and its accept/reject outcome,
  // so this digest moves if op SELECTION moves or if any op's rng consumption
  // shifts. Off unless DBARTS_FUZZ_TRACE is set, because the fuzzer's normal
  // output is one line per run; a diff of these across a change to the op
  // vocabulary is the gate that the existing configs' streams are untouched.
  if (std::getenv("DBARTS_FUZZ_TRACE") != nullptr) {
    std::uint64_t digest = 1469598103934665603ull;
    for (const std::string& text : trace)
      for (char c : text) {
        digest ^= static_cast<unsigned char>(c);
        digest *= 1099511628211ull;
      }
    printf("trace [%s] seed %u ops %zu digest %016llx\n", spec.name, seed,
           trace.size(), static_cast<unsigned long long>(digest));
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

// The multi-forest op surface: the transactional predictor paths, the
// per-observation session in both its committing and its abandoned flavor, the
// forced refresh, setCutPoints, run, the state round trip and the warm-start
// grow-from-root initializer. setData, setResponse, setWeights and setOffset
// are refused for a multi-forest sampler at the bridge
// (test-multi-forest-seam.R). OP_PER_OBS asserts, across the whole surface,
// that finalize never returns false - the guarded set is always a subset of
// the revalidated set, so an empty leaf it did not catch cannot occur - and
// OP_SESSION_ABANDON asserts that a dropped session leaves the sampler
// unchanged.
static const unsigned fuzzMultiForestMask =
  (1u << OP_SET_PREDICTOR) | (1u << OP_UPDATE_COLUMNS) | (1u << OP_PER_OBS) |
  (1u << OP_SESSION_ABANDON) | (1u << OP_SET_CUTS) | (1u << OP_RUN) |
  (1u << OP_STATE) | (1u << OP_GROW);

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
  AmplitudeSpec bcf;
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
  fuzzDrive(s, spec, arena, opRng, seed, numOps, xb, zb);

  for (size_t c = 0; c < spec.numChains; ++c) ext_rng_destroy(rngs[c]);
  ext_rng_destroy(opRng);
}

// A K-forest multinomial (softmax) sampler (docs/design/multinomial.md): the
// other shape the widened revalidation lifts, and the one with the largest
// j-splitting tree count. Its response is the borrowed count matrix, which
// OP_SET_COUNTS swaps under the running sampler and OP_SET_CATEGORY_OFFSET
// shifts. spec.family is NOT inert on this runner: fuzzDrive reads it to
// decide whether to demand a positive sigma after a run, and logistic is what
// correctly suppresses that demand on a softmax sampler.
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
                  fuzzAllSingleForestOps};
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
                    (1u << OP_SET_WEIGHTS) | (1u << OP_RUN) | (1u << OP_STATE) |
                    (1u << OP_GROW)};
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

// The structural gate on the snapshot: table-driven, one row per captured
// family, each row perturbing exactly that family and requiring the
// comparison to go false. A family with no row is a family the snapshot does
// not cover - which is the failure mode the rebuild exists to close, and which
// no amount of green fuzzing detects.
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
    AmplitudeSpec spec;
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
      {"the saved-tree draw count",
       [](FuzzSnap& t) { ++t.state.recordedDraws; }},
      {"sigma", [](FuzzSnap& t) { t.state.chains[0].sigma += 1.0; }},
      {"the BCF glue",
       [](FuzzSnap& t) { t.state.chains[0].amplitudes[0] += 1.0; }},
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

// ---------------------------------------------------------------------------
// The two gates on the per-observation session's PRUNED cache. The fuzzer
// above proves the sampler stays self-consistent whatever the session
// decided; these prove the session decided what an independent count says it
// must, and that the trees pruning skipped came through untouched to the bit.

// The three widened shapes, each owning the buffers it lends the sampler for
// its life. All instantiate Sampler<ConstantGaussianLeaf>; only the coupling
// differs, which is the point - one session implementation serves each. The
// heteroscedastic build is not multi-forest (numForests == 1) but carries the
// variance forest, the third ensemble class the session and the revalidation
// walk, so it rides the same fixture.
struct MultiForestFixture {
  std::vector<double> x, y, z;
  std::vector<int> counts, trials;
  std::vector<size_t> moderators;
  std::vector<ext_rng*> rngs;
  std::unique_ptr<Sampler<ConstantGaussianLeaf>> sampler;
  size_t n = 0, p = 0;

  ~MultiForestFixture() {
    sampler.reset();
    for (ext_rng* r : rngs) ext_rng_destroy(r);
  }

  void makeRngs(size_t numChains, std::uint32_t seed) {
    rngs.resize(numChains);
    for (size_t c = 0; c < numChains; ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], seed + static_cast<std::uint32_t>(c));
    }
  }

  // tau reads columns 0 and 2 alone, so columns 1 and 3 prune the treatment
  // forest out of the session entirely and the pruning is not vacuous
  void buildBCF(size_t n_, size_t p_, size_t numChains, std::uint32_t seed) {
    n = n_;
    p = p_;
    x.resize(n * p);
    for (double& v : x) v = runif01();
    y.resize(n);
    z.resize(n);
    for (size_t i = 0; i < n; ++i) {
      z[i] = runif01() < 0.5 ? 0.0 : 1.0;
      y[i] = 4.0 * (x[i] - 0.5) + 2.0 * x[i + n] +
             z[i] * (1.0 + 2.0 * x[i + 2 * n]) + 0.2 * (runif01() - 0.5);
    }
    moderators = {0, 2};
    SamplerOptions options;
    options.numChains = numChains;
    AmplitudeSpec spec;
    spec.mu.numTrees = 10;
    spec.tau.numTrees = 8;
    spec.tau.columns = moderators.data();
    spec.tau.numColumns = moderators.size();
    spec.z = z.data();
    makeRngs(numChains, seed);
    sampler = std::make_unique<Sampler<ConstantGaussianLeaf>>(
      x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0, fuzzRawScale,
      options, spec, rngs.data());
    Results empty;
    sampler->run(40, 0, empty);
  }

  void buildMultinomial(size_t n_, size_t p_, size_t numChains,
                        std::uint32_t seed) {
    n = n_;
    p = p_;
    const size_t K = 3;
    x.resize(n * p);
    for (double& v : x) v = runif01();
    counts.assign(n * K, 0);
    trials.assign(n, 1);
    for (size_t i = 0; i < n; ++i) {
      size_t k = runif01() < 0.3 + 0.5 * x[i] ? 0 : (runif01() < 0.5 ? 1 : 2);
      counts[k * n + i] = 1;
    }
    SamplerOptions options;
    options.numChains = numChains;
    MultinomialSpec mn;
    mn.numCategories = K;
    mn.counts = counts.data();
    mn.trials = trials.data();
    mn.forest.numTrees = 8;
    makeRngs(numChains, seed);
    sampler = std::make_unique<Sampler<ConstantGaussianLeaf>>(
      x.data(), n, p, options, mn, rngs.data());
    Results empty;
    sampler->run(40, 0, empty);
  }

  // one mean forest plus the variance forest: scale signal on column 1 and
  // mean signal on column 0, so both ensembles split and the pruning has
  // something to skip in either. restrictVariance confines the VARIANCE forest
  // to columns 0 and 1, so columns 2 and 3 are reachable by no variance tree at
  // any tree count or depth - the shape in which a transaction leaves the whole
  // variance forest untouched.
  void buildHeteroscedastic(size_t n_, size_t p_, size_t numChains,
                            std::uint32_t seed, bool restrictVariance = false) {
    n = n_;
    p = p_;
    x.resize(n * p);
    for (double& v : x) v = runif01();
    y.resize(n);
    for (size_t i = 0; i < n; ++i) {
      double u1 = runif01(), u2 = runif01();
      double e = std::sqrt(-2.0 * std::log(u1)) *
                 std::cos(6.283185307179586 * u2);
      y[i] = 3.0 * x[i] + (x[i + n] < 0.5 ? 0.25 : 1.5) * e;
    }
    SamplerOptions options;
    options.numChains = numChains;
    options.numTrees = 10;
    options.numVarianceTrees = 8;
    moderators = {0, 1};
    if (restrictVariance) {
      options.varianceForestColumns = moderators.data();
      options.numVarianceForestColumns = moderators.size();
    }
    makeRngs(numChains, seed);
    sampler = std::make_unique<Sampler<ConstantGaussianLeaf>>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, fuzzRawScale, options, rngs.data());
    Results empty;
    sampler->run(40, 0, empty);
  }
};

// An independent oracle: its own per-tree leaf counts over EVERY tree of
// EVERY forest of EVERY chain - and, on a heteroscedastic sampler, every
// variance tree - pruned or not, maintained independently of the session. The
// engine's cached set is a strict subset of this one, so a pruning bug that
// dropped a tree able to veto shows up as a decision the oracle refuses and
// the engine accepts.
struct MaskOracle {
  std::vector<const Tree*> trees;
  std::vector<std::vector<int32_t>> obsLeaf;            // tree, observation
  std::vector<std::vector<std::uint32_t>> leafCounts;   // tree, arena id

  void build(const Sampler<ConstantGaussianLeaf>& s) {
    size_t n = s.numObservations();
    trees.clear();
    for (size_t c = 0; c < s.numChains(); ++c) {
      const auto& ch(s.chain(c));
      for (size_t f = 0; f < ch.numForests(); ++f)
        for (size_t t = 0; t < ch.numTreesInForest(f); ++t)
          trees.push_back(&ch.treeInForest(f, t));
      if (!ch.hasVarianceForest()) continue;
      for (size_t j = 0; j < ch.numVarianceTrees(); ++j)
        trees.push_back(&ch.varianceTree(j));
    }
    obsLeaf.assign(trees.size(), std::vector<int32_t>(n, invalidNode));
    leafCounts.resize(trees.size());
    for (size_t k = 0; k < trees.size(); ++k) {
      leafCounts[k].assign(trees[k]->nodes.size(), 0);
      for (size_t i = 0; i < n; ++i) {
        int32_t leaf = trees[k]->findBottomNodeForObservation(s.data(), i);
        obsLeaf[k][i] = leaf;
        ++leafCounts[k][static_cast<size_t>(leaf)];
      }
    }
  }

  // the per-sampler conjunction, computed from the oracle's own counts: the
  // row installs iff it leaves every leaf of every tree occupied
  bool wouldRemainValid(const ColumnStore& data, size_t i, size_t column,
                        xint_t newCode, std::vector<int32_t>& newLeaf) const {
    bool valid = true;
    for (size_t k = 0; k < trees.size(); ++k) {
      newLeaf[k] = trees[k]->findBottomNodeForObservation(
        data, i, static_cast<int32_t>(column), newCode);
      if (newLeaf[k] != obsLeaf[k][i] &&
          leafCounts[k][static_cast<size_t>(obsLeaf[k][i])] == 1)
        valid = false;
    }
    return valid;
  }

  void commit(size_t i, const std::vector<int32_t>& newLeaf) {
    for (size_t k = 0; k < trees.size(); ++k) {
      if (newLeaf[k] == obsLeaf[k][i]) continue;
      --leafCounts[k][static_cast<size_t>(obsLeaf[k][i])];
      ++leafCounts[k][static_cast<size_t>(newLeaf[k])];
      obsLeaf[k][i] = newLeaf[k];
    }
  }
};

// One session, driven at a CALLER-CHOSEN scan order so nothing here rides the
// engine's drawn permutation. Returns the number of decisions made and adds
// their verdicts to installs/declines.
static size_t maskOracleSession(Sampler<ConstantGaussianLeaf>& s, ext_rng* r,
                                size_t column, const double* newColumn,
                                const char* label, std::vector<char>& installed,
                                size_t& installs, size_t& declines,
                                size_t& mismatches) {
  size_t n = s.numObservations();
  MaskOracle oracle;
  oracle.build(s);
  // the same quantizer the session descends on, so the comparison is about
  // the DECISION and not about how a value became a code
  std::vector<xint_t> newCodes(n);
  for (size_t i = 0; i < n; ++i)
    newCodes[i] = s.data().codeFor(column, newColumn[i]);

  std::vector<size_t> order(n);
  for (size_t i = 0; i < n; ++i) order[i] = i;
  for (size_t i = n; i > 1; --i) std::swap(order[i - 1], order[fuzzInt(r, i)]);

  std::unique_ptr<PredictorUpdateSession> session =
    s.beginPredictorUpdate(newColumn, column);
  std::vector<int32_t> newLeaf(oracle.trees.size(), invalidNode);
  installed.assign(n, 0);
  for (size_t k = 0; k < n; ++k) {
    size_t i = order[k];
    bool expected =
      oracle.wouldRemainValid(s.data(), i, column, newCodes[i], newLeaf);
    bool actual = session->observationWouldRemainValid(i);
    if (actual != expected) ++mismatches;
    if (expected) ++installs; else ++declines;
    if (actual) {
      session->commitObservation(i);
      oracle.commit(i, newLeaf);
      installed[i] = 1;
    }
  }
  // the guarded set is a subset of the revalidated set on this surface too,
  // so the finalize cannot report an empty leaf
  check(session->finalize(), label);
  return n;
}

// Candidate replacement columns. flavor 0 jitters, which almost every row can
// take; 1 collapses onto two values of the grid and 2 onto one, which empty
// leaves and force declines - so the oracle is compared on both verdicts and
// not only on the easy one.
static void maskCandidateColumn(ext_rng* r, const double* current, size_t n,
                                int flavor, std::vector<double>& out) {
  out.resize(n);
  double level = 0.2 + 0.6 * fuzzUnif(r);
  for (size_t i = 0; i < n; ++i)
    out[i] = flavor == 0   ? current[i] + 0.02 * (fuzzUnif(r) - 0.5)
             : flavor == 1 ? (fuzzUnif(r) < 0.5 ? 0.25 : 0.75)
                           : level;
}

static void testPerObservationMaskExactness() {
  ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(r, 7311u);
  const size_t n = 250, p = 4;
  size_t decisions = 0, installs = 0, declines = 0, mismatches = 0;
  size_t varianceDecisions = 0;

  for (int shape = 0; shape < 3; ++shape) {
    rngState = 0x51ED270B0F71C8ull + static_cast<uint64_t>(shape);
    MultiForestFixture fixture;
    if (shape == 0) fixture.buildBCF(n, p, 2, 4200u);
    else if (shape == 1) fixture.buildMultinomial(n, p, 1, 4300u);
    else fixture.buildHeteroscedastic(n, p, 1, 4600u);
    Sampler<ConstantGaussianLeaf>& s(*fixture.sampler);
    // the heteroscedastic shape carries the 1e5-decision floor on its own, the
    // BCF and multinomial shapes split it between them
    size_t numSessions = shape == 2 ? 420 : 220;
    // the driver keeps the current column itself; the engine holds no matrix
    std::vector<double> current(fixture.x);
    std::vector<double> candidate;
    std::vector<char> installed;
    for (size_t k = 0; k < numSessions; ++k) {
      size_t column = fuzzInt(r, p);
      int flavor = k % 8 == 7 ? 2 : (k % 4 == 3 ? 1 : 0);
      maskCandidateColumn(r, current.data() + column * n, n, flavor,
                          candidate);
      size_t made = maskOracleSession(
        s, r, column, candidate.data(),
        shape == 0 ? "the BCF session's finalize holds"
        : shape == 1 ? "the multinomial session's finalize holds"
                     : "the heteroscedastic session's finalize holds",
        installed, installs, declines, mismatches);
      decisions += made;
      if (shape == 2) varianceDecisions += made;
      // the engine keeps no predictor matrix, so the driver tracks the
      // installed cells itself, as the R layer does
      for (size_t i = 0; i < n; ++i)
        if (installed[i]) current[column * n + i] = candidate[i];
    }
  }

  check(mismatches == 0, "every install decision matches the oracle");
  check(decisions >= 100000, "at least 1e5 install decisions");
  check(varianceDecisions >= 100000,
        "at least 1e5 of them under a variance forest");
  check(installs > 0 && declines > 0, "both verdicts occur");
  ext_rng_destroy(r);
  printf("ok: per-observation mask exactness (%zu decisions, %zu under a "
         "variance forest, %zu installed, %zu declined)\n", decisions,
         varianceDecisions, installs, declines);
}

// One tree's untouched-tree families: everything an untouched tree must carry
// through a transaction unchanged.
struct F6Tree {
  std::vector<FuzzNode> nodes;   // split rules and index ranges
  std::vector<index_t> indices;  // the raw index buffer, within-leaf order too
  std::vector<double> fits;      // the tree's own fit slab
  bool operator==(const F6Tree&) const = default;
};

static bool f6FlatEqual(const std::vector<FlatNode>& a,
                        const std::vector<FlatNode>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i].variable != b[i].variable ||
        a[i].numMaskWords != b[i].numMaskWords || a[i].mask != b[i].mask ||
        a[i].flags != b[i].flags)
      return false;  // .mask is the union's raw word: a leaf value compares bitwise
  return true;
}

struct F6Capture {
  std::vector<std::vector<std::vector<F6Tree>>> trees;    // chain, forest, tree
  std::vector<std::vector<std::vector<double>>> totals;   // chain, forest
  std::vector<std::vector<F6Tree>> varianceTrees;         // chain, tree
  std::vector<std::vector<double>> combinedVariance;      // chain
  SamplerStateData state;                                 // recovered leaves
};

static F6Capture f6Capture(Sampler<ConstantGaussianLeaf>& s) {
  F6Capture out;
  size_t n = s.numObservations();
  s.getState(out.state);
  out.trees.resize(s.numChains());
  out.totals.resize(s.numChains());
  out.varianceTrees.resize(s.numChains());
  out.combinedVariance.resize(s.numChains());
  std::vector<double> slab;
  for (size_t c = 0; c < s.numChains(); ++c) {
    const auto& ch(s.chain(c));
    if (ch.hasVarianceForest()) {
      // the variance twin of the per-forest slab: the tree's own factors
      // h_j(x_i), and the product s^2(x) they compose
      const double* factors = ch.varianceFactorsForTesting();
      size_t m = ch.numVarianceTrees();
      out.varianceTrees[c].resize(m);
      for (size_t j = 0; j < m; ++j) {
        const Tree& tree(ch.varianceTree(j));
        F6Tree& e(out.varianceTrees[c][j]);
        e.nodes.resize(tree.nodes.size());
        for (size_t i = 0; i < tree.nodes.size(); ++i) {
          const Node& nd(tree.nodes[i]);
          e.nodes[i] = FuzzNode{nd.leftChild, nd.rule.variableIndex,
                                nd.rule.bits, nd.begin, nd.end};
        }
        e.indices.assign(tree.indices, tree.indices + n);
        e.fits.assign(factors + j * n, factors + (j + 1) * n);
      }
      const double* combined = ch.varianceFits();
      out.combinedVariance[c].assign(combined, combined + n);
    }
    for (size_t f = 0; f < ch.numForests(); ++f) {
      size_t numTrees = ch.numTreesInForest(f);
      slab.assign(n * numTrees, 0.0);
      ch.forestTreeFits(f, slab.data());
      std::vector<F6Tree> forestTrees(numTrees);
      for (size_t t = 0; t < numTrees; ++t) {
        const Tree& tree(ch.treeInForest(f, t));
        F6Tree& e(forestTrees[t]);
        e.nodes.resize(tree.nodes.size());
        for (size_t i = 0; i < tree.nodes.size(); ++i) {
          const Node& nd(tree.nodes[i]);
          e.nodes[i] = FuzzNode{nd.leftChild, nd.rule.variableIndex,
                                nd.rule.bits, nd.begin, nd.end};
        }
        e.indices.assign(tree.indices, tree.indices + n);
        e.fits.assign(slab.begin() + static_cast<long>(t * n),
                      slab.begin() + static_cast<long>((t + 1) * n));
      }
      out.trees[c].push_back(std::move(forestTrees));
      out.totals[c].push_back(ch.totalFitsInForest(f));
    }
  }
  return out;
}

// On the shipped (pruned) build: after an accepted transaction, a tree with
// no split on a touched column is bitwise where it was - rules, index
// ranges, the index buffer itself, its fit slab and its persisted leaf
// parameters - and if no tree of a forest splits on the column, that forest's
// totalFits keeps every contribution bitwise, where an unpruned rebuild would
// round trip it through a subtract and an add.
//
// Asserted for f >= 1 and for the variance forest, whose factor slab and
// combined variance are the scale analogues of a fit slab and a totalFits.
// Forest 0 is NEVER pruned, by rule: its subtract-then-add round trip IS the
// recorded single-forest arithmetic, and (T - f) + f differs from T at the last
// ulp on a few percent of rows, so pruning it would move an equivalence
// baseline. The asymmetry is deliberate and the count this test prints is the
// measurement of it - not a bug.
static void testUntouchedTreeExactness() {
  size_t skipped = 0, forestZeroTotalMoves = 0, varianceSkipped = 0;
  size_t varianceForestsSkipped = 0;
  for (int shape = 0; shape < 4; ++shape) {
    rngState = 0x6C1F3A55D0027Bull + static_cast<uint64_t>(shape);
    const size_t n = 240, p = 4;
    MultiForestFixture fixture;
    if (shape == 0) fixture.buildBCF(n, p, 2, 4400u);
    else if (shape == 1) fixture.buildMultinomial(n, p, 1, 4500u);
    else if (shape == 2) fixture.buildHeteroscedastic(n, p, 1, 4700u);
    // an unrestricted variance forest of any useful size splits somewhere on
    // every column, so the aggregate claim below needs the confined shape to
    // have a case at all
    else fixture.buildHeteroscedastic(n, p, 2, 4750u, true);
    Sampler<ConstantGaussianLeaf>& s(*fixture.sampler);

    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, 7411u + static_cast<std::uint32_t>(shape));
    std::vector<double> current(fixture.x);
    std::vector<std::uint32_t> uses(p);

    for (size_t column = 0; column < p; ++column) {
      // a jitter small enough to be accepted, so the transaction really runs
      // its rebuild rather than rolling back
      std::vector<double> candidate(n);
      for (size_t i = 0; i < n; ++i)
        candidate[i] = current[column * n + i] + 0.004 * (fuzzUnif(r) - 0.5);
      F6Capture before = f6Capture(s);
      PredictorUpdateResult res =
        s.updatePredictor(candidate.data(), &column, 1, false, false);
      if (res != PredictorUpdateResult::accepted) continue;
      for (size_t i = 0; i < n; ++i) current[column * n + i] = candidate[i];
      F6Capture after = f6Capture(s);

      for (size_t c = 0; c < s.numChains(); ++c) {
        const auto& ch(s.chain(c));
        for (size_t f = 1; f < ch.numForests(); ++f) {
          bool forestUntouched = true;
          for (size_t t = 0; t < ch.numTreesInForest(f); ++t) {
            std::memset(uses.data(), 0, p * sizeof(std::uint32_t));
            ch.treeInForest(f, t).countVariableUses(uses.data());
            if (uses[column] != 0) { forestUntouched = false; continue; }
            ++skipped;
            check(before.trees[c][f][t] == after.trees[c][f][t],
                  "an untouched tree's rules, ranges, indices and fits "
                  "are bitwise unchanged");
            check(f6FlatEqual(before.state.chains[c].forests[f].trees[t],
                              after.state.chains[c].forests[f].trees[t]),
                  "an untouched tree's leaf parameters are bitwise "
                  "unchanged");
          }
          if (forestUntouched)
            check(before.totals[c][f] == after.totals[c][f],
                  "a forest no touched column reaches keeps its totalFits "
                  "bitwise");
        }
        // the same for the variance forest, which prunes on the same predicate
        // with no forest-0 exemption
        if (ch.hasVarianceForest()) {
          bool varianceUntouched = true;
          for (size_t j = 0; j < ch.numVarianceTrees(); ++j) {
            std::memset(uses.data(), 0, p * sizeof(std::uint32_t));
            ch.varianceTree(j).countVariableUses(uses.data());
            if (uses[column] != 0) { varianceUntouched = false; continue; }
            ++varianceSkipped;
            check(before.varianceTrees[c][j] == after.varianceTrees[c][j],
                  "an untouched variance tree's rules, ranges, indices and "
                  "factors are bitwise unchanged");
            check(f6FlatEqual(before.state.chains[c].varianceTrees[j],
                              after.state.chains[c].varianceTrees[j]),
                  "an untouched variance tree's leaf factors are bitwise "
                  "unchanged");
          }
          if (varianceUntouched) {
            ++varianceForestsSkipped;
            check(before.combinedVariance[c] == after.combinedVariance[c],
                  "a variance forest no touched column reaches keeps "
                  "s^2(x) bitwise");
          }
        }
        // the measured asymmetry, reported rather than asserted
        for (size_t i = 0; i < n; ++i)
          if (before.totals[c][0][i] != after.totals[c][0][i])
            ++forestZeroTotalMoves;
      }
    }
    ext_rng_destroy(r);
  }
  check(skipped > 0, "the pruning actually skipped trees");
  check(varianceSkipped > 0, "the pruning skipped variance trees too");
  check(varianceForestsSkipped > 0,
        "the pruning skipped a whole variance forest");
  printf("ok: untouched-tree exactness (%zu pruned trees checked, %zu of them "
         "variance, %zu whole variance forests; forest 0's unpruned round trip "
         "moved totalFits on %zu rows)\n", skipped, varianceSkipped,
         varianceForestsSkipped, forestZeroTotalMoves);
}

// ---------------------------------------------------------------------------
// The two gates on the variance forest's half of the transaction: the
// rollback identity and the recovery ordering.

// Node-indexed leaf factors of variance tree j, read off the LIVE partition:
// the same quantity the validate phase recovers, computed independently here
// so this test compares against a value the engine did not hand us.
static void f10NodeFactors(const Tree& tree, const double* factor,
                           int32_t nodeIndex, std::vector<double>& out) {
  const Node& node(tree.at(nodeIndex));
  if (node.isBottom()) {
    if (node.numObservations() > 0)
      out[static_cast<size_t>(nodeIndex)] = factor[tree.indices[node.begin]];
    return;
  }
  f10NodeFactors(tree, factor, node.leftChild, out);
  f10NodeFactors(tree, factor, node.leftChild + 1, out);
}

// After an ACCEPTED transaction, every observation's factor is the
// PRE-transaction factor of the leaf it now lands in. That is the ordering
// claim - the recovery reads each leaf's current members before any partition
// change - stated over a quantity a caller can see. Recovering after the
// repartition instead reads the new partition's members out of the old leaves'
// slots, which moves a factor wherever the partition moved. It does not cover
// the missing-direction drop, which routes nothing and so is order-free by
// construction (tree.hpp, dropStaleMissingDirections).
static void testVarianceRecoveryOrdering() {
  rngState = 0x2D77C41B5E0093ull;
  const size_t n = 240, p = 4;
  MultiForestFixture fixture;
  fixture.buildHeteroscedastic(n, p, 1, 4800u);
  Sampler<ConstantGaussianLeaf>& s(*fixture.sampler);

  ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(r, 7511u);
  std::vector<double> current(fixture.x);
  size_t accepted = 0, moved = 0;

  for (size_t pass = 0; pass < 3; ++pass) {
    for (size_t column = 0; column < p; ++column) {
      std::vector<double> candidate(n);
      for (size_t i = 0; i < n; ++i)
        candidate[i] = current[column * n + i] + 0.01 * (fuzzUnif(r) - 0.5);

      const auto& ch(s.chain(0));
      size_t m = ch.numVarianceTrees();
      std::vector<std::vector<double>> nodeFactors(m);
      std::vector<std::vector<int32_t>> ownerBefore(m);
      const double* factorsBefore = ch.varianceFactorsForTesting();
      for (size_t j = 0; j < m; ++j) {
        const Tree& tree(ch.varianceTree(j));
        nodeFactors[j].assign(tree.nodes.size(), 0.0);
        f10NodeFactors(tree, factorsBefore + j * n, 0, nodeFactors[j]);
        ownerBefore[j].assign(n, invalidNode);
        fuzzFillOwner(tree, 0, ownerBefore[j]);
      }

      if (s.updatePredictor(candidate.data(), &column, 1, false, false) !=
          PredictorUpdateResult::accepted)
        continue;
      ++accepted;
      for (size_t i = 0; i < n; ++i) current[column * n + i] = candidate[i];

      const double* factorsAfter = ch.varianceFactorsForTesting();
      std::vector<int32_t> ownerAfter(n, invalidNode);
      for (size_t j = 0; j < m; ++j) {
        fuzzFillOwner(ch.varianceTree(j), 0, ownerAfter);
        for (size_t i = 0; i < n; ++i) {
          if (ownerAfter[i] != ownerBefore[j][i]) ++moved;
          check(factorsAfter[j * n + i] ==
                  nodeFactors[j][static_cast<size_t>(ownerAfter[i])],
                "a rebuilt factor is its new leaf's pre-transaction "
                "factor");
        }
      }
    }
  }
  check(accepted > 0, "at least one transaction was accepted");
  check(moved > 0, "the accepted transactions moved variance partitions");
  ext_rng_destroy(r);
  printf("ok: variance recovery ordering (%zu accepted transactions, %zu "
         "observations changed variance leaves)\n", accepted, moved);
}

// A ROLLED BACK transaction leaves a heteroscedastic sampler bitwise where it
// was, under the widened snapshot - which covers each variance tree's
// partition (by leaf assignment), the per-tree factor slab and s^2(x). The
// validate phase re-routes the variance trees before the veto fires, so
// repartitionTrees has to put them back; without its variance arm s^2(x) stays
// routed by the rejected proposal.
static void testVarianceRollback() {
  rngState = 0x71B9E0C33A4D15ull;
  const size_t n = 240, p = 4;
  MultiForestFixture fixture;
  fixture.buildHeteroscedastic(n, p, 1, 4900u);
  Sampler<ConstantGaussianLeaf>& s(*fixture.sampler);
  using Snap = FuzzSnapshot<Sampler<ConstantGaussianLeaf>>;

  size_t rejected = 0, varianceRerouted = 0;
  std::vector<std::uint32_t> uses(p);
  for (size_t column = 0; column < p; ++column) {
    // two levels of the existing grid: every tree splitting on the column
    // empties a leaf, so the whole transaction reverts
    std::vector<double> candidate(n);
    for (size_t i = 0; i < n; ++i) candidate[i] = i % 2 == 0 ? 0.25 : 0.75;
    // a rollback only has variance partitions to restore where some variance
    // tree splits on the column - the arm the repartition's variance step is
    // load-bearing on
    const auto& ch(s.chain(0));
    size_t varianceSplitters = 0;
    for (size_t j = 0; j < ch.numVarianceTrees(); ++j) {
      std::memset(uses.data(), 0, p * sizeof(std::uint32_t));
      ch.varianceTree(j).countVariableUses(uses.data());
      if (uses[column] != 0) ++varianceSplitters;
    }
    Snap before = fuzzCapture(s);
    PredictorUpdateResult res =
      s.updatePredictor(candidate.data(), &column, 1, false, false);
    if (res == PredictorUpdateResult::accepted) continue;
    ++rejected;
    if (varianceSplitters > 0) ++varianceRerouted;
    check(fuzzSnapshotsEqual(before, fuzzCapture(s)),
          "a rolled-back transaction leaves the variance forest bitwise");
  }
  // the whole-matrix flavor too, whose null column list prunes nothing
  {
    std::vector<double> candidate(n * p);
    for (size_t i = 0; i < n * p; ++i) candidate[i] = i % 2 == 0 ? 0.25 : 0.75;
    Snap before = fuzzCapture(s);
    if (s.setPredictor(candidate.data(), false, false) !=
        PredictorUpdateResult::accepted) {
      ++rejected;
      check(fuzzSnapshotsEqual(before, fuzzCapture(s)),
            "a rolled-back whole-matrix swap leaves the variance forest "
            "bitwise");
    }
  }
  check(rejected > 0, "at least one transaction rolled back");
  check(varianceRerouted > 0,
        "at least one rollback had variance partitions to restore");
  printf("ok: variance rollback identity (%zu rolled-back transactions, %zu "
         "with variance partitions to restore)\n", rejected, varianceRerouted);
}

static void testMutationFuzzer(int numSeeds) {
  ColSpec ord;
  ColSpec cat4{ColumnType::categorical, 4, 0.0};
  ColSpec ordMiss{ColumnType::ordinal, 0, 0.25};
  ColSpec cat3Miss{ColumnType::categorical, 3, 0.2};
  std::vector<ConfigSpec> configs = {
    {"gaussian", ResponseFamily::gaussian, {ord, ord, ord}, 1,
     fuzzAllSingleForestOps},
    {"categorical", ResponseFamily::gaussian, {cat4, ord, ord}, 1,
     fuzzAllSingleForestOps},
    {"missing", ResponseFamily::gaussian, {ordMiss, cat3Miss, ord}, 1,
     fuzzAllSingleForestOps},
    {"multichain", ResponseFamily::gaussian, {ord, ord, ord}, 2,
     fuzzAllSingleForestOps},
    {"probit", ResponseFamily::probit, {ord, ord, ord}, 1,
     (1u << OP_SET_PREDICTOR) | (1u << OP_UPDATE_COLUMNS) | (1u << OP_PER_OBS) |
       (1u << OP_SESSION_ABANDON) | (1u << OP_SET_DATA) | (1u << OP_SET_CUTS) |
       (1u << OP_SET_TEST) | (1u << OP_RUN) | (1u << OP_STATE) |
       (1u << OP_GROW)},
    {"logistic", ResponseFamily::logistic, {ord, ord, ord}, 1,
     (1u << OP_SET_PREDICTOR) | (1u << OP_UPDATE_COLUMNS) | (1u << OP_PER_OBS) |
       (1u << OP_SESSION_ABANDON) | (1u << OP_SET_DATA) | (1u << OP_SET_CUTS) |
       (1u << OP_SET_TEST) | (1u << OP_RUN) | (1u << OP_STATE) |
       (1u << OP_GROW)},
    // setCuts joined the safe surface once setCutPoints re-routed the variance
    // forest through forceRefreshTrees, and setData once applyNewData resized
    // its seven n-sized allocations and re-routed it - the op the fuzzer runs
    // at a fresh count every time. The transactional predictor paths and the
    // per-observation session joined once the two-phase revalidation and the
    // session's cell guard reached the variance forest, and OP_GROW once
    // grow-from-root's own formMeanWeights branch (chain.hpp) was confirmed
    // to scan against the variance forest's precisions rather than unit
    // ones, so this config now runs every op but setSigma, which stays out
    // because it is an engine no-op under a variance forest (the forest IS
    // the residual variance).
    {"heteroscedastic", ResponseFamily::gaussian, {ord, ord, ord}, 1,
     (1u << OP_SET_RESPONSE) | (1u << OP_SET_WEIGHTS) | (1u << OP_SET_OFFSET) |
       (1u << OP_SET_TEST) | (1u << OP_RUN) | (1u << OP_STATE) |
       (1u << OP_SET_CUTS) | (1u << OP_SET_DATA) |
       (1u << OP_SET_PREDICTOR) | (1u << OP_UPDATE_COLUMNS) |
       (1u << OP_PER_OBS) | (1u << OP_SESSION_ABANDON) | (1u << OP_GROW),
     8},
  };
  // the multi-forest shapes, on their own runners: one continuous and one
  // categorical design each, so the transactional paths meet both column types
  ConfigSpec bcf{"bcf", ResponseFamily::gaussian, {ord, ord, ord}, 1,
                 fuzzMultiForestMask};
  ConfigSpec bcfCat{"bcf-categorical", ResponseFamily::gaussian,
                    {ord, cat4, ord}, 2, fuzzMultiForestMask};
  // the multinomial rows carry the response-side ops on top of that mask;
  // the two-chain row is what drives Sampler::setCounts' and
  // setCategoryOffset's per-chain fan-out and I1's chain stride
  ConfigSpec multinomial{"multinomial", ResponseFamily::logistic,
                         {ord, ord, ord}, 1,
                         fuzzMultiForestMask | fuzzCountsOps};
  ConfigSpec multinomialMiss{"multinomial-missing", ResponseFamily::logistic,
                             {ordMiss, cat3Miss, ord}, 1,
                             fuzzMultiForestMask | fuzzCountsOps};
  ConfigSpec multinomialChains{"multinomial-multichain",
                               ResponseFamily::logistic, {ord, ord, ord}, 2,
                               fuzzMultiForestMask | fuzzCountsOps};

  const int numOps = 40;
  for (int sd = 0; sd < numSeeds; ++sd) {
    std::uint32_t seed = 1u + static_cast<std::uint32_t>(sd);
    for (const ConfigSpec& spec : configs)
      fuzzRunConstant(spec, seed, numOps);
    fuzzRunBCF(bcf, seed, numOps);
    fuzzRunBCF(bcfCat, seed, numOps);
    fuzzRunMultinomial(multinomial, seed, numOps);
    fuzzRunMultinomial(multinomialMiss, seed, numOps);
    fuzzRunMultinomial(multinomialChains, seed, numOps);
    fuzzRunLinear(seed, numOps);
    fuzzRunSparse(seed, numOps);
  }
  printf("ok: mutation fuzzer (%d seeds)\n", numSeeds);
}

void runFuzzTests(int numSeeds) {
  testSnapshotCoversEveryFamily();
  testPerObservationMaskExactness();
  testUntouchedTreeExactness();
  testVarianceRecoveryOrdering();
  testVarianceRollback();
  testMutationFuzzer(numSeeds);
}
