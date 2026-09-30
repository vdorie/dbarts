// The monotone "leaf" prior's birth/death ratio Z_T0 / Z_T*, counted on the
// side of the finer tree T* (docs/plans/monotone-exact-birth-death.md,
// Counting: algorithm), against the direct count of the merged component C0.
//
// Input on stdin, one move per line, as benchmarks/R/monotone-order-size.R
// writes it (pairs mode):
//   pair <id> <split variable> <m> <c1> <c2> <m masks> | <n0> <n0 masks>
// The first order holds the components of T* that contain the two children c1
// and c2 (c1 the lower one), the second the component C0 of T0 that holds
// their merged leaf. Each element is two 64-bit hex words (high, low), its
// predecessor mask (bit j of element k set when j < k), as monotone_count
// reads them.
//
// Per move it prints: id, "one" or "two" (components of T* holding the pair),
// their sizes and down-set counts, log theta and seconds on T*'s side; then
// C0's down-set count (by a branching recursion that needs no layers) and,
// when that is within the cap, the direct log theta = log e(C0) - log e(U) and
// its seconds. theta is the probability that c2 immediately follows c1 in a
// uniform linear extension of T*'s components, and Z_T0 / Z_T* = m theta.
//
// With a second argument K > 0 it also draws K uniform linear extensions of
// T*'s components by Huber's bounding-chain coupling from the past (Huber
// 2006, Discrete Math. 306: 420-428), one component at a time and riffled,
// and prints the share where c2 immediately follows c1 and the milliseconds
// per draw: the coin of the count-free alternative in the plan.
//
// Usage: ./monotone_ratio [cap on direct down-sets, default 1e7] [K, default 0]
// Build: make monotone_ratio (needs no package build).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <utility>
#include <vector>

typedef unsigned __int128 u128;
typedef std::vector<u128> Order;  // predecessor masks, transitively closed

struct Hash {
  size_t operator()(u128 k) const {
    uint64_t a = (uint64_t) k, b = (uint64_t) (k >> 64);
    return a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL + (a << 6));
  }
};

static double now() {
  return std::chrono::duration<double>(
           std::chrono::steady_clock::now().time_since_epoch())
    .count();
}
static inline bool has(u128 m, int j) { return (m >> j) & 1; }
static inline u128 one(int j) { return (u128) 1 << j; }

// ---- counts --------------------------------------------------------------

// Layered down-set DP, two layers: log e(P), or NAN past the cap.
static double logExtensions(const Order& pred, double cap, double* downSets) {
  int n = (int) pred.size();
  std::unordered_map<u128, double, Hash> cur, next;
  cur[0] = 1.0;
  double total = 1.0;
  for (int k = 0; k < n; ++k) {
    next.clear();
    next.reserve(cur.size() * 2);
    for (const auto& kv : cur)
      for (int x = 0; x < n; ++x)
        if (!has(kv.first, x) && !(pred[x] & ~kv.first))
          next[kv.first | one(x)] += kv.second;
    total += (double) next.size();
    if (total > cap) {
      *downSets = total;
      return NAN;
    }
    std::swap(cur, next);
  }
  *downSets = total;
  return std::log(cur.begin()->second);
}

// Every layer, with forward counts f(D) (orders of D) and backward counts g(D)
// (linear extensions of the rest): P(x at position i) sums f(D) g(D + x) over
// the down-sets D of size i that x can extend.
static std::vector<double> positionLaw(const Order& pred, int x, double* ds) {
  int n = (int) pred.size();
  std::vector<std::unordered_map<u128, std::pair<double, double>, Hash>> layer(
    n + 1);
  layer[0][0] = {1.0, 0.0};
  *ds = 1.0;
  for (int k = 0; k < n; ++k) {
    layer[k + 1].reserve(layer[k].size() * 2);
    for (const auto& kv : layer[k])
      for (int y = 0; y < n; ++y)
        if (!has(kv.first, y) && !(pred[y] & ~kv.first))
          layer[k + 1][kv.first | one(y)].first += kv.second.first;
    *ds += (double) layer[k + 1].size();
  }
  layer[n].begin()->second.second = 1.0;
  for (int k = n - 1; k >= 0; --k)
    for (auto& kv : layer[k]) {
      double g = 0.0;
      for (int y = 0; y < n; ++y)
        if (!has(kv.first, y) && !(pred[y] & ~kv.first))
          g += layer[k + 1].find(kv.first | one(y))->second.second;
      kv.second.second = g;
    }
  double e = layer[0].begin()->second.second;
  std::vector<double> law(n, 0.0);
  for (int k = 0; k < n; ++k)
    for (const auto& kv : layer[k])
      if (!has(kv.first, x) && !(pred[x] & ~kv.first))
        law[k] += kv.second.first *
                  layer[k + 1].find(kv.first | one(x))->second.second / e;
  return law;
}

// Down-sets by branching on one element (in or out), memoized on what is
// left: cheap on these orders even where the layered DP is not.
static std::unordered_map<u128, double, Hash> memo;
static Order above, below;
static double downSetsOf(u128 set) {
  if (!set) return 1.0;
  auto hit = memo.find(set);
  if (hit != memo.end()) return hit->second;
  int x = 0;
  while (!has(set, x)) ++x;
  double r = downSetsOf(set & ~one(x) & ~above[x]) +
             downSetsOf(set & ~one(x) & ~below[x]);
  memo[set] = r;
  return r;
}
static double countDownSets(const Order& pred) {
  int n = (int) pred.size();
  below = pred;
  above.assign(n, 0);
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      if (has(pred[k], j)) above[j] |= one(k);
  memo.clear();
  return downSetsOf(n == 128 ? ~(u128) 0 : one(n) - 1);
}

static Order restrict(const Order& pred, const std::vector<int>& ix) {
  Order out(ix.size(), 0);
  for (size_t a = 0; a < ix.size(); ++a)
    for (size_t b = 0; b < ix.size(); ++b)
      if (has(pred[ix[b]], ix[a])) out[b] |= one((int) a);
  return out;
}

// c2 merged into c1 (the union of their relations, closed): T0's order
static Order mergePair(const Order& pred, int c1, int c2) {
  int n = (int) pred.size();
  std::vector<int> keep;
  for (int i = 0; i < n; ++i)
    if (i != c2) keep.push_back(i);
  Order u = pred;  // fold c2 into c1 before restricting
  for (int k = 0; k < n; ++k) {
    if (has(u[k], c2)) u[k] |= one(c1);
  }
  u[c1] |= u[c2];
  u[c1] &= ~(one(c1) | one(c2));
  Order out = restrict(u, keep);
  for (bool changed = true; changed;) {
    changed = false;
    for (size_t b = 0; b < out.size(); ++b) {
      u128 m = out[b];
      for (size_t a = 0; a < out.size(); ++a)
        if (has(out[b], (int) a)) m |= out[a];
      if (m != out[b]) {
        out[b] = m;
        changed = true;
      }
    }
  }
  return out;
}

static std::vector<std::vector<int>> components(const Order& pred) {
  int n = (int) pred.size();
  std::vector<int> label(n, -1);
  std::vector<std::vector<int>> out;
  for (int s = 0; s < n; ++s) {
    if (label[s] >= 0) continue;
    std::vector<int> stack{s}, comp;
    label[s] = (int) out.size();
    while (!stack.empty()) {
      int x = stack.back();
      stack.pop_back();
      comp.push_back(x);
      for (int y = 0; y < n; ++y)
        if (label[y] < 0 && (has(pred[y], x) || has(pred[x], y))) {
          label[y] = (int) out.size();
          stack.push_back(y);
        }
    }
    std::sort(comp.begin(), comp.end());
    out.push_back(comp);
  }
  return out;
}

static double choose[129][129];

// ---- Huber (2006): bounding-chain CFTP for a uniform linear extension -----

struct Rng {
  uint64_t s[2];
  explicit Rng(uint64_t seed) {
    s[0] = seed * 0x9E3779B97F4A7C15ULL + 1;
    s[1] = (seed ^ 0xDEADBEEFCAFEULL) * 0xBF58476D1CE4E5B9ULL + 7;
    for (int i = 0; i < 10; ++i) next();
  }
  uint64_t next() {
    uint64_t s1 = s[0];
    const uint64_t s0 = s[1];
    s[0] = s0;
    s1 ^= s1 << 23;
    s[1] = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return s[1] + s0;
  }
};

// pred must have the identity as a linear extension; returns the elements in
// order. The Karzanov-Khachiyan chain (hold 1/2, else a uniform adjacent
// transposition when allowed) runs coupled to the bounding chain R, where
// X^-1(a) <= R(a) for the elements in play; the coupling flips the coin when
// X(i) is the in-play element with R = i + 1 (Huber's Theorem 3), so it is
// non-Markovian and each level's forward pass reruns the bounding chain.
static std::vector<int> huberDraw(const Order& pred, Rng& master) {
  int n = (int) pred.size();
  if (n == 1) return {0};
  std::vector<int> R(n), who(n), X(n);
  auto run = [&](uint64_t seed, long len, bool advanceX) {
    Rng g(seed);
    std::fill(R.begin(), R.end(), n - 1);
    std::fill(who.begin(), who.end(), -1);
    who[n - 1] = 0;
    int inPlay = 1;
    for (long t = 0; t < len; ++t) {
      uint64_t u = g.next();
      int i = (int) ((u >> 1) % (uint64_t) (n - 1));
      int c = (int) (u & 1);
      int a = who[i], b = who[i + 1];
      if (advanceX) {
        int cx = (b >= 0 && X[i] == b) ? 1 - c : c;
        if (cx && !has(pred[X[i + 1]], X[i])) std::swap(X[i], X[i + 1]);
      }
      if (c) {
        if (a >= 0 && b >= 0) {
          if (!has(pred[b], a)) {
            R[a] = i + 1;
            R[b] = i;
            who[i] = b;
            who[i + 1] = a;
          }
        } else if (a >= 0) {
          R[a] = i + 1;
          who[i] = -1;
          who[i + 1] = a;
        } else if (b >= 0) {
          R[b] = i;
          who[i + 1] = -1;
          who[i] = b;
        }
      }
      if (inPlay < n && who[n - 1] < 0) {
        R[inPlay] = n - 1;
        who[n - 1] = inPlay++;
      }
    }
    return inPlay == n;  // every element in play: R is a bijection
  };
  std::vector<uint64_t> seeds;
  long base = std::max(16L, (long) n * n);
  int level = 0;
  for (;; ++level) {
    seeds.push_back(master.next());
    if (run(seeds.back(), base << level, false)) break;
  }
  for (int a = 0; a < n; ++a) X[R[a]] = a;
  for (int l = level - 1; l >= 0; --l) run(seeds[l], base << l, true);
  return X;
}

// a uniform linear extension of pred (any labels): relabel topologically
static std::vector<int> drawExtension(const Order& pred, Rng& rng) {
  int n = (int) pred.size();
  std::vector<int> ord(n);
  for (int i = 0; i < n; ++i) ord[i] = i;
  auto size = [&](int k) {
    return __builtin_popcountll((uint64_t) pred[k]) +
           __builtin_popcountll((uint64_t) (pred[k] >> 64));
  };
  std::stable_sort(ord.begin(), ord.end(),
                   [&](int a, int b) { return size(a) < size(b); });
  std::vector<int> label(n);
  for (int i = 0; i < n; ++i) label[ord[i]] = i;
  Order p(n, 0);
  for (int b = 0; b < n; ++b)
    for (int a = 0; a < n; ++a)
      if (has(pred[b], a)) p[label[b]] |= one(label[a]);
  std::vector<int> x = huberDraw(p, rng);
  for (int& v : x) v = ord[v];
  return x;
}

// ---- main ------------------------------------------------------------------

static Order readOrder(int n, char*& s) {
  Order pred(n);
  for (int k = 0; k < n; ++k) {
    unsigned long long hi = strtoull(s, &s, 16), lo = strtoull(s, &s, 16);
    pred[k] = ((u128) hi << 64) | lo;
  }
  return pred;
}

int main(int argc, char** argv) {
  double cap = argc > 1 ? atof(argv[1]) : 1e7;
  int draws = argc > 2 ? atoi(argv[2]) : 0;
  for (int n = 0; n <= 128; ++n) {
    choose[n][0] = choose[n][n] = 1.0;
    for (int k = 1; k < n; ++k)
      choose[n][k] = choose[n - 1][k - 1] + choose[n - 1][k];
  }
  Rng rng(20260929);
  static char line[1 << 17];
  double worstRel = 0.0, zSum = 0.0, z2Sum = 0.0;
  int checked = 0, zN = 0;
  std::printf(
    "id split n1 n2 D1 D2 logTheta sec D0 directLogTheta directSec%s\n",
    draws ? " drawShare msPerDraw" : "");
  while (fgets(line, sizeof line, stdin)) {
    char id[160];
    int var, m, c1, c2, off;
    char* s = line;
    if (sscanf(s, "pair %159s %d %d %d %d%n", id, &var, &m, &c1, &c2, &off) !=
        5)
      continue;
    s += off;
    Order U = readOrder(m, s);
    while (*s == ' ' || *s == '|') ++s;
    int n0 = (int) strtol(s, &s, 10);
    Order C0 = readOrder(n0, s);

    std::vector<std::vector<int>> comps = components(U);
    int k1 = -1, k2 = -1;
    for (int k = 0; k < (int) comps.size(); ++k)
      for (int x : comps[k]) {
        if (x == c1) k1 = k;
        if (x == c2) k2 = k;
      }
    auto at = [](const std::vector<int>& v, int x) {
      return (int) (std::find(v.begin(), v.end(), x) - v.begin());
    };
    double t0 = now(), logTheta, D1, D2 = 0.0;
    int n1 = (int) comps[k1].size(), n2 = 0;
    if (k1 == k2) {
      // one component C*: theta = e(C0) / e(C*), and D(C0) <= D(C*)
      Order Cs = restrict(U, comps[k1]);
      Order merged = mergePair(Cs, at(comps[k1], c1), at(comps[k1], c2));
      logTheta = logExtensions(merged, INFINITY, &D2) -
                 logExtensions(Cs, INFINITY, &D1);
    } else {
      // two components: position laws of c1 in C1 and c2 in C2, and the
      // share of riffles that put the one just before the other
      n2 = (int) comps[k2].size();
      std::vector<double> p1 = positionLaw(restrict(U, comps[k1]),
                                           at(comps[k1], c1), &D1);
      std::vector<double> p2 = positionLaw(restrict(U, comps[k2]),
                                           at(comps[k2], c2), &D2);
      double theta = 0.0;
      for (int i = 1; i <= n1; ++i)
        for (int j = 1; j <= n2; ++j)
          theta += p1[i - 1] * p2[j - 1] * choose[i + j - 2][i - 1] *
                   choose[n1 + n2 - i - j][n1 - i];
      logTheta = std::log(theta / choose[n1 + n2][n1]);
    }
    double sec = now() - t0;

    // the plan's direct count: e(C0) and e of the union of T*'s components
    double D0 = countDownSets(C0), directLog = NAN, directSec = NAN;
    if (D0 <= cap) {
      double t1 = now(), d;
      double logEU = std::lgamma(m + 1.0);
      std::vector<int> held{k1};
      if (k2 != k1) held.push_back(k2);
      for (int k : held) {
        Order Ck = restrict(U, comps[k]);
        logEU +=
          logExtensions(Ck, INFINITY, &d) - std::lgamma(comps[k].size() + 1.0);
      }
      directLog = logExtensions(C0, INFINITY, &d) - logEU;
      directSec = now() - t1;
      worstRel = std::max(worstRel, std::fabs(std::expm1(logTheta - directLog)));
      ++checked;
    }
    std::printf("%s %s %d %d %.4g %.4g %.12g %.3g %.4g %.12g %.3g", id,
                k1 == k2 ? "one" : "two", n1, n2, D1, D2, logTheta, sec, D0,
                directLog, directSec);

    if (draws) {
      int hits = 0;
      double t2 = now();
      for (int d = 0; d < draws; ++d) {
        std::vector<int> seq;
        std::vector<int> a = drawExtension(restrict(U, comps[k1]), rng);
        for (int& v : a) v = comps[k1][v];
        if (k1 == k2) {
          seq = a;
        } else {
          std::vector<int> b = drawExtension(restrict(U, comps[k2]), rng);
          for (int& v : b) v = comps[k2][v];
          size_t ia = 0, ib = 0;
          while (ia < a.size() || ib < b.size()) {
            uint64_t left = a.size() - ia, total = left + b.size() - ib;
            if (rng.next() % total < left)
              seq.push_back(a[ia++]);
            else
              seq.push_back(b[ib++]);
          }
        }
        int p1 = at(seq, c1), p2 = at(seq, c2);
        hits += p2 == p1 + 1;
      }
      double share = (double) hits / draws, theta = std::exp(logTheta);
      std::printf(" %.4g %.3g", share, 1e3 * (now() - t2) / draws);
      if (theta > 0.0 && theta < 1.0) {
        double z = (share - theta) / std::sqrt(theta * (1 - theta) / draws);
        zSum += z;
        z2Sum += z * z;
        ++zN;
      }
    }
    std::printf("\n");
    std::fflush(stdout);
  }
  std::printf("# %d moves checked against the direct count: worst relative "
              "difference in theta %.3g\n",
              checked, worstRel);
  if (zN)
    std::printf("# draws: %d moves, mean z %.3f, mean z^2 %.3f\n", zN,
                zSum / zN, z2Sum / zN);
  return 0;
}
