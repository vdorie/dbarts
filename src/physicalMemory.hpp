#ifndef DBARTS_PHYSICAL_MEMORY_HPP
#define DBARTS_PHYSICAL_MEMORY_HPP

#include <cstdint>

namespace dbarts {
  /// Returns the physical memory installed in the machine, in bytes, or 0
  /// where the operating system does not report it.
  ///
  /// This is what the machine has, not what this process may use: neither the
  /// memory free at the call nor a limit set by a container or a job scheduler
  /// is read here, so under such a limit the figure is larger than what the
  /// process can allocate. The R caller reads a Linux cgroup limit itself.
  std::uint64_t physicalMemory() noexcept;
}

#endif // DBARTS_PHYSICAL_MEMORY_HPP
