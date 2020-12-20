#include "config.hpp"
#include "guessNumCores.hpp"

#include <cstddef> // size_t

#include <external/io.h>

using std::size_t;
using std::uint32_t;

// windows adapted from
//   http://scriptionary.com/2010/09/01/counting-processor-cores-and-threads/
//   http://msdn.microsoft.com/en-us/library/ms683194%28v=VS.85%29.aspx
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>

#  ifndef _W64
#    include "glpi.h"
#  endif

typedef BOOL (WINAPI* glpiFunction)(PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, PDWORD);

namespace dbarts {
  void guessNumCores(uint32_t* numPhysicalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhysicalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    glpiFunction getLogicalProcessorInformation = (glpiFunction) GetProcAddress(GetModuleHandle(TEXT("kernel32")), "GetLogicalProcessorInformation");
    if (getLogicalProcessorInformation == NULL) {
      ext_issueWarning("GetLogicalProcessorInformation is not supported on this OS.");
      return;
    }
    
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = NULL;
    DWORD bufferLength = 0;
    
    // retrieve appropriate buffer length
    getLogicalProcessorInformation(buffer, &bufferLength);
    if (GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
      ext_issueWarning("Unable to determine buffer length for GetLogicalProcessorInformation: %d", GetLastError());
      return;
    }
    
    buffer = new SYSTEM_LOGICAL_PROCESSOR_INFORMATION[bufferLength];
    if (getLogicalProcessorInformation(buffer, &bufferLength) == false) {
      delete [] buffer;
      ext_issueWarning("Unable to allocate buffer length for GetLogicalProcessorInformation");
      return;
    }
    
    size_t numElements = bufferLength / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
    
    for (size_t i = 0; i < numElements; ++i) {
      if (buffer[i].Relationship != RelationProcessorCore) continue;
      
      *numPhysicalProcessorsPtr += 1;
      
      if (buffer[i].ProcessorCore.Flags == 1) {
        ULONG mask = 0x1;
        size_t numBitsInLong = 8 * sizeof(ULONG);
        for (size_t j = 0; j < numBitsInLong; ++j) {
          if ((buffer[i].ProcessorMask & mask) != 0) *numLogicalProcessorsPtr += 1;
          mask <<= 1;
        }
      }
    }
    
    delete [] buffer;
  }
}

#elif (defined(__APPLE__) || defined(__FreeBSD__) || defined (__NetBSD__))
#  include <sys/types.h>
#  include <sys/sysctl.h>

namespace dbarts {
  void guessNumCores(uint32_t* numPhysicalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhysicalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    int query[2];
    size_t queryLength = 2;
    
    size_t bufferLength = sizeof(uint32_t);
    
    // these have apparently been around since 2003 so I'm not sure when/if they could fail
    if (sysctlnametomib("hw.physicalcpu", query, &queryLength) != -1) {
      sysctl(query, static_cast<u_int>(queryLength), numPhysicalProcessorsPtr, &bufferLength, NULL, 0);
    }
    
    queryLength = 2;
    bufferLength = sizeof(uint32_t);
    
    if (sysctlnametomib("hw.logicalcpu", query, &queryLength) != -1) {
      sysctl(query, static_cast<u_int>(queryLength), numLogicalProcessorsPtr, &bufferLength, NULL, 0);
    }
    
    bufferLength = sizeof(uint32_t);
    
    if (*numLogicalProcessorsPtr == 0) {
      query[0] = CTL_HW;
#  ifdef HW_AVAILCPU
      query[1] = HW_AVAILCPU;
      if (sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0) == -1 || *numLogicalProcessorsPtr == 0) {
#  endif
        query[1] = HW_NCPU;
        sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0);
      }
#  ifdef HW_AVAILCPU
    }
#  endif
  }
}

// the only way to reliably determine the number of physical cores
// is to parse /proc/cpuinfo; sysconf(_SC_NPROCESSORS_CONF) returns
// hyper as well
#elif defined(__linux__)
#  include <fcntl.h>    // open
#  include <unistd.h>   // close, lseek
#  include <sys/mman.h> // mmap
#  include <cstring>   // srtncmp, strerror_r
#  include <errno.h>
#  include <sys/stat.h> // stat

#  include <misc/alloca.h>

#  include <map>
#  include <vector>
#  include <utility>

#  include <cstdio>
#  include <cstdlib>

namespace {
  typedef std::map<uint32_t, uint32_t> CoreMap;
  typedef std::pair<uint32_t, uint32_t> CorePair;
  struct Processor {
    uint32_t physicalId;
    CoreMap coreIdThreadMap;
  };
  typedef std::map<uint32_t, Processor*> CPUMap;
  typedef std::pair<uint32_t, Processor*> CPUPair;
  
  enum LineType {
    PROCESSOR = 0,
    PHYSICAL_ID,
    CORE_ID,
    NA
  };
  
#  define INVALID_ID ((uint32_t) -1)
  
  void parseRecord(CPUMap& processorIdMap, uint32_t& currentPhysicalId, uint32_t& currentCoreId, uint32_t& currentProcessorId, uint32_t newId)
  {
    Processor* currentProcessor;
    
    // new records are introduced by processor id lines
    if (currentProcessorId != newId && currentProcessorId != INVALID_ID) {
      
      if (currentPhysicalId == INVALID_ID) {
        // only should happen when there is a single processor with one core
        // and no hyperthreading
        currentProcessor = new Processor;
        currentProcessor->physicalId = 0;
        currentProcessor->coreIdThreadMap.insert(CorePair(0, 1));
        
        processorIdMap.insert(CPUPair(currentProcessorId, currentProcessor));
      } else {
        // the following line should never happen
        if (currentCoreId == INVALID_ID) currentCoreId = 0;
        
        CPUMap::iterator processorLookup = processorIdMap.find(currentPhysicalId);
        if (processorLookup == processorIdMap.end()) {
          currentProcessor = new Processor;
          currentProcessor->physicalId = currentPhysicalId;
          currentProcessor->coreIdThreadMap.insert(CorePair(currentCoreId, 1));
          processorIdMap.insert(CPUPair(currentPhysicalId, currentProcessor));
        } else {
          currentProcessor = processorLookup->second;
          CoreMap::iterator coreLookup = currentProcessor->coreIdThreadMap.find(currentCoreId);
          if (coreLookup == currentProcessor->coreIdThreadMap.end()) {
            currentProcessor->coreIdThreadMap.insert(CorePair(currentCoreId, 1));
          } else {
            coreLookup->second += 1;
          }
        }
      }
      
      
      currentProcessorId = newId;
      currentPhysicalId = INVALID_ID;
      currentCoreId = INVALID_ID;
    } else {
      // first line, usually
      currentProcessorId = newId;
    }
  }

#  define ERROR_BUFFER_LENGTH static_cast<size_t>(1024)
  bool parseProcCPUInfo(std::vector < Processor* >& result)
  {
    char errorBuffer[ERROR_BUFFER_LENGTH];
    
    size_t bufferLength = sysconf(_SC_PAGESIZE);
    char* cpuInfo = new char[bufferLength];
    
    int fd = open("/proc/cpuinfo", O_RDONLY);
    if (fd == -1) {
#ifdef _GNU_SOURCE
      char* errorMessage = strerror_r(errno, errorBuffer, ERROR_BUFFER_LENGTH);
      ext_issueWarning("unable to open /proc/cpuinfo: %s (%d)\n", errorMessage, errno);
#else
      int ignored = strerror_r(errno, errorBuffer, ERROR_BUFFER_LENGTH);
      ext_issueWarning("unable to open /proc/cpuinfo: %s (%d)\n", errorBuffer, errno);
#endif
      return false;
    }
    
    // Ideally, this would simply read as it goes but I wrote it on a system w/o
    // procfs and didn't know that I can't mmap the sucker.
    
    ssize_t numBytesRead = 0;
    size_t fileLength = 0, numBytesInBuffer = 0;
    bool firstReallocation = true;
    while (true) {
      numBytesRead = read(fd, cpuInfo + fileLength, bufferLength - numBytesInBuffer);
      if (numBytesRead <= 0) break;
      fileLength += numBytesRead;
      numBytesInBuffer += numBytesRead;
      
      if (numBytesInBuffer == bufferLength) {
        char* temp = new char[2 * fileLength];
        std::memcpy(temp, (const char*) cpuInfo, fileLength * sizeof(char));
        
        if (!firstReallocation) bufferLength *= 2;
        else firstReallocation = false;
        
        numBytesInBuffer = 0;
        delete [] cpuInfo;
        cpuInfo = temp;
      }
    }
    close(fd);
    if (numBytesRead == -1) {
#ifdef _GNU_SOURCE
      char* errorMessage = strerror_r(errno, errorBuffer, ERROR_BUFFER_LENGTH);
      ext_issueWarning("unable to open /proc/cpuinfo: %s (%d)\n", errorMessage, errno);
#else
      int ignored = strerror_r(errno, errorBuffer, ERROR_BUFFER_LENGTH);
      ext_issueWarning("unable to open /proc/cpuinfo: %s (%d)\n", errorBuffer, errno);
#endif
      return false;
    }
    
    
    
    off_t offset = 0;
    
    uint32_t currentProcessorId = INVALID_ID;
    uint32_t currentPhysicalId  = INVALID_ID;
    uint32_t currentCoreId      = INVALID_ID;
    CPUMap processorIdMap;
    LineType lineType;
    
    while (offset < (off_t) fileLength) {
      lineType = NA;
      
      if ((cpuInfo[offset] == 'p' || cpuInfo[offset] == 'P') && fileLength - offset >= 9 &&
          std::strncmp(cpuInfo + offset + 1, "rocessor", 8) == 0) {
        lineType = PROCESSOR;
        offset += 9;
      } else if ((cpuInfo[offset] == 'p' || cpuInfo[offset] == 'P') && fileLength - offset >= 11 &&
                 std::strncmp(cpuInfo + offset + 1, "hysical id", 10) == 0) {
        lineType = PHYSICAL_ID;
        offset += 11;
      } else if ((cpuInfo[offset] == 'c' || cpuInfo[offset] == 'C') && fileLength - offset >= 7 &&
                 std::strncmp(cpuInfo + offset + 1, "ore id", 6) == 0) {
        lineType = CORE_ID;
        offset += 7;
      }
      
      if (lineType == NA) {
        while (offset < (off_t) fileLength && cpuInfo[offset] != '\n') ++offset;
        ++offset;
        if (offset >= (off_t) fileLength) break;
        continue;
      }
      
      while (offset < (off_t) fileLength && cpuInfo[offset] != ':' && cpuInfo[offset] != '\n') ++offset;
      if (cpuInfo[offset] == '\n') {
        ++offset; continue;
      }
      if (offset == (off_t) fileLength) break;
      ++offset;
      
      while (offset < (off_t) fileLength && cpuInfo[offset] < '0' && cpuInfo[offset] > '9' && cpuInfo[offset] != '\n') ++offset;
      if (cpuInfo[offset] == '\n') {
        ++offset; continue;
      }
      if (offset == (off_t) fileLength) break;
      
      off_t endOfNumber = offset + 1;
      while (endOfNumber < (off_t) fileLength && cpuInfo[endOfNumber] >= '0' && cpuInfo[endOfNumber] <= '9') ++endOfNumber;
      
      char* buffer = misc_stackAllocate(endOfNumber - offset + 1, char);
      std::memcpy(buffer, (const char*) cpuInfo + offset, endOfNumber - offset);
      buffer[endOfNumber - offset] = '\0';
      
      long parsedInt = std::strtol(buffer, NULL, 10);
      misc_stackFree(buffer);
      
      offset = endOfNumber;
      
      
      switch (lineType) {
        case PROCESSOR:
          parseRecord(processorIdMap, currentPhysicalId, currentCoreId, currentProcessorId, parsedInt);
          break;
        case PHYSICAL_ID:
          currentPhysicalId = parsedInt;
          break;
        case CORE_ID:
          currentCoreId = parsedInt;
          break;
        default:
          break;
      }
      
      while (offset < (off_t) fileLength && cpuInfo[offset] != '\n') ++offset;
      ++offset;
      if (offset >= (off_t) fileLength) break;
    }
    
    parseRecord(processorIdMap, currentPhysicalId, currentCoreId, currentProcessorId, INVALID_ID);
    
    for (CPUMap::iterator it = processorIdMap.begin(); it != processorIdMap.end(); ++it) {
      result.push_back(it->second);
    }
    
    delete [] cpuInfo;
    
    return true;
  }
}

namespace dbarts {
  void guessNumCores(uint32_t* numPhyiscalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhyiscalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    std::vector<Processor*> cpuInfo;
    if (parseProcCPUInfo(cpuInfo) == true) {
      for (size_t i = 0; i < cpuInfo.size(); ++i) {
        Processor* processor = cpuInfo[i];
        *numPhyiscalProcessorsPtr += processor->coreIdThreadMap.size();
        for (CoreMap::iterator it = processor->coreIdThreadMap.begin(); it != processor->coreIdThreadMap.end(); ++it) {
          *numLogicalProcessorsPtr += it->second;
        }
      }
    }
#  if defined(_SC_NPROCESSORS_ONLN)
    else {
      ext_printf("  cpuinfo falling back\n");
      *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_ONLN);
      if (*numLogicalProcessorsPtr < 1) *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_CONF);
    }
#  endif
    
    for (size_t i = 0; i < cpuInfo.size(); ++i) {
      Processor* processor = cpuInfo[i];
      delete processor;
    }
  }
}

#else // freeball it

#  include <unistd.h>
#  ifdef _SC_NPROCESSORS_ONLN
#    define USE_SYSCONF
#  elif defined(HAVE_SYS_SYSCTL_H)
#    include <sys/param.h>
#    include <sys/sysctl.h>
#    define USE_SYSCTL
#  endif


namespace dbarts {
  void guessNumCores(uint32_t* numPhysicalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhysicalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    
#  if defined(USE_SYSCONF)
    *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_ONLN);
    if (*numLogicalProcessorsPtr < 1) {
#    ifdef _SC_NPROCESSORS_CONF
      *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_CONF);
      if (*numLogicalProcessorsPtr < 1) *numLogicalProcessorsPtr = -1;
#    else
      *numLogicalProcessorsPtr = -1;
#    endif
    }
    
#  elif defined(USE_SYSCTL)
    int query[2];
    size_t bufferLength = sizeof(uint32_t);
    
    query[0] = CTL_HW;
#    ifdef HW_AVAILCPU
    query[1] = HW_AVAILCPU;
    if (sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0) == -1 ||  numCores < 1) {
#    endif
      key[1] = HW_NCPU;
      sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0);
#    ifdef HW_AVAILCPU
    }
#    endif
#  endif
  }
}

#endif

