#ifndef DBARTS_GUESS_NUM_CORES
#define DBARTS_GUESS_NUM_CORES

#include <cstdint>

namespace dbarts {
  void guessNumCores(std::uint32_t* numPhysicalProcessors, std::uint32_t* numLogicalProcessors);
}

#endif

