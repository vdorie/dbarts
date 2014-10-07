#ifndef DBARTS_BINARY_IO_HPP
#define DBARTS_BINARY_IO_HPP

struct ext_binaryIO;

namespace dbarts {
  struct Control;
  struct Data;
  struct Model;
  struct Scratch;
  struct State;
  
  bool writeControl(const Control& control, ext_binaryIO* bio);
  bool readControl(Control& control, ext_binaryIO* bio);

  bool writeData(const Data& data, ext_binaryIO* bio);
  bool readData(Data& data, ext_binaryIO* bio);
  
  bool writeModel(const Model& model, ext_binaryIO* bio);
  bool readModel(Model& model, ext_binaryIO* bio);
  
  bool writeState(const State& state, ext_binaryIO* bio, const Control& control, const Data& data);
  bool readState(State& state, ext_binaryIO* bio, const Control& control, const Data& data);
}

#endif // DBARTS_BINARY_IO_HPP
