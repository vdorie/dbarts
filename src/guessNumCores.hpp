#ifndef DBARTS_GUESS_NUM_CORES
#define DBARTS_GUESS_NUM_CORES

#include <dbarts/cstdint.hpp>

namespace dbarts {
  void guessNumCores(std::uint32_t* numPhyiscalProcessors, std::uint32_t* numLogicalProcessors);
}

#endif

