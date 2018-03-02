#ifndef DBARTS_FUNCTIONS_HPP
#define DBARTS_FUNCTIONS_HPP

#include <cstddef>
#include <dbarts/cstdint.hpp>

#include <dbarts/types.hpp>

namespace dbarts {
  struct BARTFit;
  struct State;
  struct Node;
  struct Tree;
  
  void updateVariablesAvailable(const BARTFit& fit, Node& node, std::int32_t variableIndex);
  double metropolisJumpForTree(const BARTFit& fit, std::size_t chainNum, Tree& tree, const double* y, double sigma,
                               bool* stepTaken, StepType* stepType);
  
  // for int "ind", flips the bits in d such that d is the powers of 2 to get back to ind
  void setBinaryRepresentation(std::uint32_t length, std::uint32_t ind, bool* d);
  
  std::size_t countTrueValues(bool* v, std::size_t length);
  
  int32_t findIndexOfIthPositiveValue(bool* values, std::size_t numValues, std::size_t i);
  void setCategoryReachability(const BARTFit& fit, const Node& node, std::int32_t variableIndex, bool* categoriesCanReachNode);
  void setSplitInterval(const BARTFit& fit, const Node& startNode, std::int32_t variableIndex, std::int32_t* leftIndex, std::int32_t* rightIndex);
}

#endif
