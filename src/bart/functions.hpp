#ifndef BART_FUNCTIONS_HPP
#define BART_FUNCTIONS_HPP

#include <cstddef>
#include "bart/cstdint"

#include <bart/types.hpp>

using std::size_t;
using std::int32_t;
using std::uint32_t;

namespace bart {
  struct BARTFit;
  struct Node;
  struct Tree;
  
  void updateVariablesAvailable(const BARTFit& fit, Node& node, int32_t variableIndex);
  double metropolisJumpForTree(const BARTFit& fit, Tree& tree, const double* y,
                               bool* stepTaken, StepType* stepType);
  
  // for int "ind", flips the bits in d such that d is the powers of 2 to get back to ind
  void setBinaryRepresentation(uint32_t length, uint32_t ind, bool* d);
  
  size_t countTrueValues(bool* v, size_t length);
  
  int32_t findIndexOfIthPositiveValue(bool* values, size_t numValues, size_t i);
  void setCategoryReachability(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* categoriesCanReachNode);
  void setSplitInterval(const BARTFit& fit, const Node& startNode, int32_t variableIndex, int32_t* leftIndex, int32_t* rightIndex);
}

#endif
