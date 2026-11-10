#include "physicalMemory.hpp"

#include <cstddef> // size_t
#include <cstdint> // uint64_t

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#elif defined(__APPLE__)
#  include <sys/types.h>
#  include <sys/sysctl.h>
#else
#  include <unistd.h>
#endif

using std::uint64_t;

namespace dbarts {

#if defined(_WIN32)

uint64_t physicalMemory() noexcept {
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (GlobalMemoryStatusEx(&status) == 0) return 0;
  return static_cast<uint64_t>(status.ullTotalPhys);
}

#elif defined(__APPLE__)

uint64_t physicalMemory() noexcept {
  uint64_t numBytes = 0;
  std::size_t length = sizeof(numBytes);
  if (sysctlbyname("hw.memsize", &numBytes, &length, nullptr, 0) != 0 ||
      length != sizeof(numBytes))
    return 0;
  return numBytes;
}

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

// Linux, and the other POSIX systems that report their page count the same way
uint64_t physicalMemory() noexcept {
  const long numPages = sysconf(_SC_PHYS_PAGES);
  const long pageSize = sysconf(_SC_PAGESIZE);
  // a failed query returns -1
  if (numPages <= 0 || pageSize <= 0) return 0;
  return static_cast<uint64_t>(numPages) * static_cast<uint64_t>(pageSize);
}

#else

uint64_t physicalMemory() noexcept { return 0; }

#endif

} // namespace dbarts
