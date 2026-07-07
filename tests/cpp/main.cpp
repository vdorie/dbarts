// Component and smoke tests for bartcore against independently coded
// reference math. Exit code 0 on success.
//
// Usage: test_bartcore [seed-count] [suite]
//   seed-count  an integer: the mutation fuzzer's seed count (default 3)
//   suite       a name prefix (data, tree, moves, model, sampler, state,
//               fuzz) selecting only that suite; both arguments may appear
//               in either order

#include "common.hpp"

namespace {

bool isNumeric(const char* arg, long* value) {
  char* end = nullptr;
  *value = std::strtol(arg, &end, 10);
  return end != arg && *end == '\0';
}

bool suiteSelected(const char* filter, const char* suite) {
  if (filter == nullptr) return true;
  size_t len = std::strlen(filter);
  return len > 0 && std::strncmp(filter, suite, len) == 0;
}

}  // namespace

int main(int argc, char** argv) {
  misc_simd_init();

  int numFuzzSeeds = 3;
  const char* filter = nullptr;
  for (int i = 1; i < argc; ++i) {
    long value;
    if (isNumeric(argv[i], &value)) numFuzzSeeds = static_cast<int>(value);
    else filter = argv[i];
  }
  if (numFuzzSeeds < 1) numFuzzSeeds = 1;

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rng == NULL || ext_rng_setSeed(rng, 42) != 0) {
    printf("FAIL: could not create rng\n");
    return 1;
  }

  if (suiteSelected(filter, "data")) runDataTests();
  if (suiteSelected(filter, "tree")) runTreeTests(rng);
  if (suiteSelected(filter, "moves")) runMovesTests(rng);
  if (suiteSelected(filter, "model")) runModelTests(rng);
  if (suiteSelected(filter, "sampler")) runSamplerTests(rng);
  if (suiteSelected(filter, "state")) runStateTests(rng);

  ext_rng_destroy(rng);

  if (suiteSelected(filter, "fuzz")) runFuzzTests(numFuzzSeeds);

  if (failures > 0) {
    printf("%d failure(s)\n", failures);
    return 1;
  }
  printf("all tests passed\n");
  return 0;
}
