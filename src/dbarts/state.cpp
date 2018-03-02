#include "config.hpp"
#include <dbarts/state.hpp>

#include <algorithm> // integer min
#include <cerrno>
#include <cstdlib>
#include <cstring>
#ifdef HAVE_STD_SNPRINTF
// snprintf in c++11, before that have to use C version
#  include <cstdio>
using std::snprintf;
#else
#  include <stdio.h>
#endif

#include <external/io.h>
#include <external/linearAlgebra.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include "functions.hpp"
#include "node.hpp"
#include "tree.hpp"

#define BASE_BUFFER_SIZE 1024
#define INT_BUFFER_SIZE 16

using std::size_t;

namespace {
  using namespace dbarts;
  
  void setNewObservationIndices(Node& newNode, size_t* indices, const Node& oldNode)
  {
    newNode.setObservationIndices(indices);
    if (!newNode.isBottom()) {
      setNewObservationIndices(*newNode.getLeftChild(), indices, *oldNode.getLeftChild());
      setNewObservationIndices(*newNode.getRightChild(), indices + oldNode.getLeftChild()->getNumObservations(), *oldNode.getRightChild());
    }
  }
}

namespace dbarts {
  State::State(const Control& control, const Data& data) {
    size_t numSamples = control.runMode == FIXED_SAMPLES ? control.defaultNumSamples : 1;
    
    size_t totalNumTrees = control.numTrees * numSamples;
    
    treeIndices = new size_t[data.numObservations * totalNumTrees];
    
    trees = static_cast<Tree*>(::operator new (totalNumTrees * sizeof(Tree)));
    
    for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum)
      new (trees + treeNum) Tree(treeIndices + treeNum * data.numObservations, data.numObservations, data.numPredictors);
    
    
    treeFits  = new double[data.numObservations * totalNumTrees];
    ext_setVectorToConstant(treeFits, data.numObservations * totalNumTrees, 0.0);
    
    sigma = new double[numSamples];
    
    rng = NULL;
  }
  
  void State::invalidate(size_t totalNumTrees) {
    delete [] sigma;
    delete [] treeFits;
    
    for (size_t treeNum = totalNumTrees; treeNum > 0; --treeNum)
      trees[treeNum - 1].~Tree();
    
    delete [] treeIndices;
    
    ::operator delete (trees);
  }
}

namespace {
  struct ResizeData {
    const Data& data;
    
    const Control& oldControl;
    const State& oldState;
    
    const Control& newControl;
    State& newState;
  };
    
  void copyTreesForSample(ResizeData& resizeData, size_t oldSampleNum, size_t newSampleNum) {
    const Data& data(resizeData.data);
    
    const Control& oldControl(resizeData.oldControl);
    const State& oldState(resizeData.oldState);
    const Control& newControl(resizeData.newControl);
    State& newState(resizeData.newState);
    
    size_t assignEnd = std::min(oldControl.numTrees, newControl.numTrees);
    
    // copy in trees that will be in the new one
    for (size_t treeNum = 0; treeNum < assignEnd; ++treeNum) {
      size_t oldTreeOffset = treeNum + oldSampleNum * oldControl.numTrees;
      size_t newTreeOffset = treeNum + newSampleNum * newControl.numTrees;
      
      // copy over children pointers, if any
      newState.trees[newTreeOffset] = oldState.trees[oldTreeOffset];
      setNewObservationIndices(newState.trees[newTreeOffset].top,
                               newState.treeIndices + newTreeOffset * data.numObservations,
                               oldState.trees[oldTreeOffset].top);
      
      if (!newState.trees[newTreeOffset].top.isBottom()) {
        newState.trees[newTreeOffset].top.getRightChild()->parent = &newState.trees[newTreeOffset].top;
        newState.trees[newTreeOffset].top.getLeftChild()->parent  = &newState.trees[newTreeOffset].top;
        
        // prevent destructor from freeing children since we just assigned them over
        oldState.trees[oldTreeOffset].top.leftChild = NULL;
      }
    }
    
    size_t oldSampleOffset = oldSampleNum * data.numObservations * oldControl.numTrees;
    size_t newSampleOffset = newSampleNum * data.numObservations * newControl.numTrees;
    
    std::memcpy(newState.treeIndices + newSampleOffset, oldState.treeIndices + oldSampleOffset, assignEnd * data.numObservations * sizeof(size_t));
    std::memcpy(newState.treeFits    + newSampleOffset, oldState.treeFits    + oldSampleOffset, assignEnd * data.numObservations * sizeof(double));
    
    // if any new trees are required, create and initialize those
    for (size_t treeNum = assignEnd; treeNum < newControl.numTrees; ++treeNum) {
      size_t treeOffset = treeNum + newSampleNum * newControl.numTrees;
      new (newState.trees + treeOffset)
        Tree(newState.treeIndices + treeOffset * data.numObservations, data.numObservations, data.numPredictors);
      ext_setVectorToConstant(newState.treeFits + treeOffset * data.numObservations, data.numObservations, 0.0);
    }
    
    // if any extra trees exist, delete them
    oldSampleOffset = oldSampleNum * oldControl.numTrees;
    
    for (size_t treeNum = oldControl.numTrees; treeNum > assignEnd; --treeNum)
      oldState.trees[treeNum - 1 + oldSampleOffset].~Tree();
  }
}

namespace dbarts {
  bool State::resize(const BARTFit& fit, const Control& newControl) {
    const Control& oldControl(fit.control);
    const Data& data(fit.data);
   
    if (oldControl.runMode == newControl.runMode && oldControl.numTrees == newControl.numTrees) return false;
        
    size_t oldNumSamples = oldControl.runMode == FIXED_SAMPLES ? fit.currentNumSamples : 1;
    size_t newNumSamples = newControl.runMode == FIXED_SAMPLES ? fit.currentNumSamples : 1;
    
    State oldState = *this;
    
    // trees first
    trees       = static_cast<Tree*>(::operator new (newControl.numTrees * newNumSamples * sizeof(Tree)));
    treeIndices = new size_t[data.numObservations * newControl.numTrees * newNumSamples];
    treeFits    = new double[data.numObservations * newControl.numTrees * newNumSamples];
    
    ResizeData resizeData = { fit.data, oldControl, oldState, newControl, *this };
    
    size_t numSamplesToCopy, oldSampleStart, newSampleStart;
    if (oldNumSamples > newNumSamples) {
      numSamplesToCopy = 1;
      oldSampleStart = oldNumSamples - 1;
      newSampleStart = 0;
      
      for (size_t sampleNum = oldNumSamples - 1; sampleNum > 0; --sampleNum) {
        for (size_t treeNum = oldControl.numTrees; treeNum > 0; --treeNum) {
          size_t treeOffset = (treeNum - 1) + (sampleNum - 1) * oldControl.numTrees;
          oldState.trees[treeOffset].~Tree();
        }
      }
    } else if (oldNumSamples < newNumSamples) {
      numSamplesToCopy = 1;
      oldSampleStart = 0;
      newSampleStart = newNumSamples - 1;
      
      for (size_t sampleNum = 0; sampleNum < newNumSamples - 1; ++sampleNum) {
        size_t sampleOffset = sampleNum * newControl.numTrees;
        for (size_t treeNum = 0; treeNum < newControl.numTrees; ++treeNum) {
          new (trees + treeNum + sampleOffset)
            Tree(treeIndices + data.numObservations * (treeNum + sampleOffset), data.numObservations, data.numPredictors);
        }
      }
    } else {
      numSamplesToCopy = newNumSamples;
      oldSampleStart = 0;
      newSampleStart = 0;
    }
    
    //ext_printf("copying trees for resize\n");
    for (size_t sampleNum = 0; sampleNum < numSamplesToCopy; ++sampleNum) {
      copyTreesForSample(resizeData, oldSampleStart + sampleNum, newSampleStart + sampleNum);
      //for (size_t treeNum = 0; treeNum < newControl.numTrees; ++treeNum)
      //  trees[treeNum + (sampleNum + newSampleStart) * newControl.numTrees].top.checkIndices(fit, trees[treeNum + (sampleNum + newSampleStart) * newControl.numTrees].top);
    }
    //ext_printf("tree copy done\n");
    
    
    ::operator delete (oldState.trees);
    delete [] oldState.treeIndices;
    delete [] oldState.treeFits;
    
    if (oldNumSamples != newNumSamples) {
      sigma = new double[newNumSamples];
      std::memcpy(sigma + newSampleStart, oldState.sigma + oldSampleStart, numSamplesToCopy * sizeof(double));
      delete [] oldState.sigma;
    }
    
    return true;
  }
  
  bool State::resize(const BARTFit& fit, size_t newNumSamples) {
    const Control& control(fit.control);
    const Data& data(fit.data);
    
    size_t oldNumSamples = fit.currentNumSamples;
    
    if (newNumSamples == oldNumSamples || control.runMode == SEQUENTIAL_UPDATES) return false;
    
    State oldState = *this;
    
    // trees first
    trees       = static_cast<Tree*>(::operator new (control.numTrees * newNumSamples * sizeof(Tree)));
    treeIndices = new size_t[data.numObservations * control.numTrees * newNumSamples];
    treeFits    = new double[data.numObservations * control.numTrees * newNumSamples];
    
    ResizeData resizeData = { data, control, oldState, control, *this };
    
    size_t numSamplesToCopy, oldSampleStart, newSampleStart;
    if (oldNumSamples > newNumSamples) {
      numSamplesToCopy = newNumSamples;
      oldSampleStart = oldNumSamples - newNumSamples;
      newSampleStart = 0;
      
      for (size_t sampleNum = oldNumSamples - newNumSamples; sampleNum > 0; --sampleNum) {
        for (size_t treeNum = control.numTrees; treeNum > 0; --treeNum) {
          size_t treeOffset = (treeNum - 1) + (sampleNum - 1) * control.numTrees;
          oldState.trees[treeOffset].~Tree();
        }
      }
    } else {
      numSamplesToCopy = oldNumSamples;
      oldSampleStart = 0;
      newSampleStart = newNumSamples - oldNumSamples;
      
      for (size_t sampleNum = 0; sampleNum < newNumSamples - oldNumSamples; ++sampleNum) {
        size_t sampleOffset = sampleNum * control.numTrees;
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          new (trees + treeNum + sampleOffset)
            Tree(treeIndices + data.numObservations * (treeNum + sampleOffset), data.numObservations, data.numPredictors);
        }
      }
    }
        
    for (size_t sampleNum = 0; sampleNum < numSamplesToCopy; ++sampleNum)
      copyTreesForSample(resizeData, oldSampleStart + sampleNum, newSampleStart + sampleNum);
    
     ::operator delete (oldState.trees);
    delete [] oldState.treeIndices;
    delete [] oldState.treeFits;
    
    sigma = new double[newNumSamples];
    std::memcpy(sigma + newSampleStart, oldState.sigma + oldSampleStart, numSamplesToCopy * sizeof(double));
    delete [] oldState.sigma;
    
    return true;
  }
}


namespace {
  using namespace dbarts;
  
  struct StringWriter {
    char* buffer;
    size_t length;
    size_t pos;
    
    void writeChar(char c) {
      buffer[pos++] = c;
      
      if (pos >= length) reallocate();
    }
    void writeString(const char* s, size_t len) {
      if (pos + len >= length) reallocate();
      std::memcpy(buffer + pos, s, len * sizeof(char));
      pos += len;
    }
    
    void writeInt(int32_t i) {
      char intBuffer[INT_BUFFER_SIZE];
      int bytesWritten = snprintf(intBuffer, INT_BUFFER_SIZE, "%d", i);
      writeString(intBuffer, static_cast<size_t>(bytesWritten));
    }
    
    /* void writeUInt(uint32_t u) {
      char intBuffer[INT_BUFFER_SIZE];
      int bytesWritten = snprintf(intBuffer, INT_BUFFER_SIZE, "%u", u);
      writeString(intBuffer, (size_t) bytesWritten);
    } */
    
    void writeNode(const Node& node) {
      if (node.isBottom()) {
        writeChar('.');
        return;
      }
      
      writeInt(node.p.rule.variableIndex);
      writeChar(' ');
      writeInt(node.p.rule.splitIndex);
      writeChar(' ');
      
      writeNode(*node.getLeftChild());
      writeNode(*node.getRightChild());
    }
    
    void reallocate() {
      char* temp = new char[length + BASE_BUFFER_SIZE];
      std::memcpy(temp, const_cast<const char*>(buffer), length * sizeof(char));
      
      length += BASE_BUFFER_SIZE;
      delete [] buffer;
      buffer = temp;
    }
  };

  using namespace dbarts;
  
  size_t readNode(Node& node, const char* treeString, size_t numPredictors) {
    if (treeString[0] == '\0') return 0;
    if (treeString[0] == '.') return 1;
    
    
    size_t pos = 0;
    
    char buffer[INT_BUFFER_SIZE];
    while (treeString[pos] != ' ' && pos < INT_BUFFER_SIZE) {
      buffer[pos] = treeString[pos];
      ++pos;
    }
    
    if (pos == INT_BUFFER_SIZE) ext_throwError("Unable to parse tree string: expected integer.");
    buffer[pos++] = '\0';
    
    
    errno = 0;
    node.p.rule.variableIndex = static_cast<int32_t>(std::strtol(buffer, NULL, 10));
    if (node.p.rule.variableIndex == 0 && errno != 0)
      ext_throwError("Unable to parse tree string: %s", std::strerror(errno));
    
    size_t bufferPos = 0;
    while (treeString[pos] != ' ' && bufferPos < INT_BUFFER_SIZE) {
      buffer[bufferPos++] = treeString[pos++];
    }
    
    if (pos == INT_BUFFER_SIZE) ext_throwError("Unable to parse tree string: expected integer.");
    buffer[bufferPos++] = '\0';
    ++pos;
    
    errno = 0;
    node.p.rule.splitIndex = static_cast<int32_t>(std::strtol(buffer, NULL, 10));
    if (node.p.rule.splitIndex == 0 && errno != 0)
      ext_throwError("Unable to parse tree string: %s", std::strerror(errno));
    
    node.leftChild  = new Node(node, numPredictors);
    node.p.rightChild = new Node(node, numPredictors);
    
    pos += readNode(*node.getLeftChild(), treeString + pos, numPredictors);
    pos += readNode(*node.getRightChild(), treeString + pos, numPredictors);
    
    return pos;
  }
}


namespace dbarts {
  const char* const* State::createTreeStrings(const BARTFit& fit) const
  {
    StringWriter writer;
    
    size_t numTrees = fit.control.numTrees * (fit.control.runMode == FIXED_SAMPLES ? fit.currentNumSamples : 1);
    
    char** result = new char*[numTrees];
    for (size_t i = 0; i < numTrees; ++i) {
      writer.buffer = new char[BASE_BUFFER_SIZE];
      writer.length = BASE_BUFFER_SIZE;
      writer.pos = 0;
      
      writer.writeNode(trees[i].top);
     
      writer.writeChar('\0');
      
      result[i] = writer.buffer;
    }
    
    return(result);
  }
  
  void State::recreateTreesFromStrings(const BARTFit& fit, const char* const* treeStrings)
  {
    size_t numTrees = fit.control.numTrees * (fit.control.runMode == FIXED_SAMPLES ? fit.currentNumSamples : 1);
    
    for (size_t i = 0; i < numTrees; ++i) {
      trees[i].top.clear();
      readNode(trees[i].top, treeStrings[i], fit.data.numPredictors);
      
      if (!trees[i].top.isBottom()) {
        updateVariablesAvailable(fit, trees[i].top, trees[i].top.p.rule.variableIndex);
      
        trees[i].top.addObservationsToChildren(fit);
      }
    }
  }
}
