// Layered down-set DP for the monotone normalizer Z_T = e(P_T) / L!
// (docs/plans/monotone-exact-birth-death.md, Decision): times the count the
// planned engine runs per component of a tree's leaf order, keeping two layers
// of a hash map from down-set bitset to the number of linear extensions of it.
//
// Input on stdin, one order per line: an id, the element count n <= 128, then
// n pairs of 64-bit hex words (high, low) giving each element's predecessor
// mask (bit j of element k set when j < k). benchmarks/R/monotone-order-size.R
// writes this format for a fit's components (masks=<file>).
//
// Output, one line per order: id n down-sets ok|over log(e) seconds, then a
// '#' line with the total seconds per down-set x element. A count stops once
// its down-sets pass the cap (first argument, default 1e8).
//
// `./monotone_count synthetic` times two worst-case shapes at about 2^24
// down-sets instead: a 25-element star (one minimum below 24 others) and a
// 60-element order with one minimum below eight chains (lengths 8,8,8,7,7,7,7,7).
//
// Build: make monotone_count (needs no package build).

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

typedef unsigned __int128 u128;

struct Hash {
  size_t operator()(u128 k) const {
    uint64_t a = (uint64_t) k, b = (uint64_t) (k >> 64);
    return a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL + (a << 6));
  }
};

static double totalSeconds = 0.0, totalUnits = 0.0;

static void count(const char* id, const std::vector<u128>& pred, double cap) {
  int n = (int) pred.size();
  auto t0 = std::chrono::steady_clock::now();
  std::unordered_map<u128, double, Hash> cur, next;
  cur[0] = 1.0;
  double downSets = 1.0;
  bool over = false;
  for (int layer = 0; layer < n && !over; ++layer) {
    next.clear();
    next.reserve(cur.size() * 2);
    for (const auto& kv : cur) {
      u128 d = kv.first;
      for (int x = 0; x < n; ++x) {
        u128 bit = (u128) 1 << x;
        if ((d & bit) || (pred[x] & ~d)) continue;
        next[d | bit] += kv.second;
      }
    }
    downSets += (double) next.size();
    over = downSets > cap;
    std::swap(cur, next);
  }
  double seconds =
    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
  double logE = over ? NAN : std::log(cur.begin()->second);
  totalSeconds += seconds;
  totalUnits += downSets * n;
  std::printf(
    "%s %d %.0f %s %.6g %.6g\n", id, n, downSets, over ? "over" : "ok", logE,
    seconds);
  std::fflush(stdout);
}

int main(int argc, char** argv) {
  bool synthetic = argc > 1 && std::strcmp(argv[1], "synthetic") == 0;
  double cap = argc > 1 && !synthetic ? std::atof(argv[1]) : 1e8;
  if (synthetic) {
    std::vector<u128> star(25, (u128) 1);
    star[0] = 0;
    count("star25", star, cap);
    std::vector<u128> chains(1, (u128) 0);
    const int lengths[] = {8, 8, 8, 7, 7, 7, 7, 7};
    for (int len : lengths) {
      u128 below = 1;  // the shared minimum
      for (int i = 0; i < len; ++i) {
        chains.push_back(below);
        below |= (u128) 1 << (chains.size() - 1);
      }
    }
    count("chains60", chains, cap);
  } else {
    char id[256];
    int n;
    while (std::scanf("%255s %d", id, &n) == 2) {
      std::vector<u128> pred(n);
      for (int k = 0; k < n; ++k) {
        unsigned long long hi, lo;
        if (std::scanf("%llx %llx", &hi, &lo) != 2) return 1;
        pred[k] = ((u128) hi << 64) | lo;
      }
      count(id, pred, cap);
    }
  }
  std::printf(
    "# %.3g ns per down-set x element\n", 1e9 * totalSeconds / totalUnits);
  return 0;
}
