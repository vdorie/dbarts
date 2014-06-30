#ifndef DBARTS_TYPES_HPP
#define DBARTS_TYPES_HPP

namespace dbarts {
  enum VariableType {
    ORDINAL, CATEGORICAL
  };
  
  enum StepType {
    BIRTH, DEATH, SWAP, CHANGE
  };
} // namespace dbarts

#endif // DBARTS_TYPES_HPP
