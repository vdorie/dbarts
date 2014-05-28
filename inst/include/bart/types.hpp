#ifndef BART_TYPES_HPP
#define BART_TYPES_HPP

namespace bart {
  enum VariableType {
    ORDINAL, CATEGORICAL
  };
  
  enum StepType {
    BIRTH, DEATH, SWAP, CHANGE
  };
} // namespace bart

#endif // BART_TYPES_HPP
