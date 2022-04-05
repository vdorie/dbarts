#ifndef DBARTS_BINARY_IO_HPP
#define DBARTS_BINARY_IO_HPP

struct misc_binaryIO;

#include <cstddef>

namespace dbarts {
  struct Control;
  struct Data;
  struct Model;
  struct State;
  
  struct Version {
    std::size_t major;
    std::size_t minor;
    std::size_t revision;
  };
  
  bool writeControl(misc_binaryIO* bio, const Control& control);
  bool readControl(misc_binaryIO* bio, Control& control, const Version& version);

  bool writeData(misc_binaryIO* bio, const Data& data);
  bool readData(misc_binaryIO* bio, Data& data);
  
  bool writeModel(misc_binaryIO* bio, const Model& model);
  bool readModel(misc_binaryIO* bio, Model& model);
  
  bool writeState(misc_binaryIO* bio, const State* state, const Control& control, const Data& data, std::size_t numSamples, std::size_t treeFitsStride);
  bool readState(misc_binaryIO* bio, State* state, const Control& control, const Data& data, std::size_t numSamples, std::size_t treeFitsStride, const Version& version);
  
  int readVersion(misc_binaryIO* bio, Version& version, char** versionStringPtr);
}

#endif // DBARTS_BINARY_IO_HPP
