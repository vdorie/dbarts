#ifndef DBARTS_BINARY_IO_HPP
#define DBARTS_BINARY_IO_HPP

struct ext_binaryIO;

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
  
  bool writeControl(ext_binaryIO* bio, const Control& control);
  bool readControl(ext_binaryIO* bio, Control& control, const Version& version);

  bool writeData(ext_binaryIO* bio, const Data& data);
  bool readData(ext_binaryIO* bio, Data& data);
  
  bool writeModel(ext_binaryIO* bio, const Model& model);
  bool readModel(ext_binaryIO* bio, Model& model);
  
  bool writeState(ext_binaryIO* bio, const State* state, const Control& control, const Data& data);
  bool readState(ext_binaryIO* bio, State* state, const Control& control, const Data& data, const Version& version);
  
  int readVersion(ext_binaryIO* bio, Version& version, char** versionStringPtr);
}

#endif // DBARTS_BINARY_IO_HPP
