#ifndef DBARTS_TYPES_HPP
#define DBARTS_TYPES_HPP

#include "config.hpp"

namespace dbarts {
  enum VariableType {
    ORDINAL, CATEGORICAL
  };
  
  enum StepType {
    BIRTH, DEATH, SWAP, CHANGE
  };
  
  typedef std::XINT_TYPE xint_t;
} // namespace dbarts

#endif // DBARTS_TYPES_HPP
