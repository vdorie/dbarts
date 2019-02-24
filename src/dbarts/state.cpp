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
#include <misc/linearAlgebra.h>

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
  State::State(const Control& control, const Data& data)
  {
    size_t totalNumTrees = control.numTrees;
    
    treeIndices = new size_t[data.numObservations * totalNumTrees];
    
    trees = static_cast<Tree*>(::operator new (totalNumTrees * sizeof(Tree)));
    for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum)
      new (trees + treeNum) Tree(treeIndices + treeNum * data.numObservations, data.numObservations, data.numPredictors);
    
    treeFits = new double[data.numObservations * totalNumTrees];
    misc_setVectorToConstant(treeFits, data.numObservations * totalNumTrees, 0.0);
    
    if (control.keepTrees) {
      totalNumTrees *= control.defaultNumSamples;
      
      savedTreeIndices = new size_t[data.numObservations * totalNumTrees];
      
      savedTrees = static_cast<Tree*>(::operator new (totalNumTrees * sizeof(Tree)));
      for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum)
        new (savedTrees + treeNum) Tree(savedTreeIndices + treeNum * data.numObservations, data.numObservations, data.numPredictors);
      
      savedTreeFits = new double[data.numObservations * totalNumTrees];
      misc_setVectorToConstant(savedTreeFits, data.numObservations * totalNumTrees, 0.0);
    } else {
      savedTreeIndices = NULL;
      savedTrees = NULL;
      savedTreeFits = NULL;
    }
        
    rng = NULL;
  }
  
  void State::invalidate(size_t numTrees, size_t numSamples) {
    delete [] savedTreeFits;
    if (savedTrees != NULL) {
      for (size_t treeNum = numTrees * numSamples; treeNum > 0; --treeNum)
        savedTrees[treeNum - 1].~Tree();
      ::operator delete (savedTrees);
    }
    delete [] savedTreeIndices;

    
    delete [] treeFits;
    for (size_t treeNum = numTrees; treeNum > 0; --treeNum)
      trees[treeNum - 1].~Tree();
    ::operator delete (trees);
    delete [] treeIndices;
  }
}

namespace {
  struct TreeData {
    std::size_t* treeIndices; // numObs x numTrees
    Tree* trees;              // numTrees
    double* treeFits;
  };
  struct ResizeData {
    const Data& data;
    
    const Control& oldControl;
    const Control& newControl;
    
    const TreeData& oldTreeData;
    TreeData& newTreeData;
  };
    
  void copyTreesForSample(ResizeData& resizeData, size_t oldSampleNum, size_t newSampleNum) {
    const Data& data(resizeData.data);
    
    const Control& oldControl(resizeData.oldControl);
    const Control& newControl(resizeData.newControl);
    
    const TreeData& oldTreeData(resizeData.oldTreeData);
    TreeData& newTreeData(resizeData.newTreeData);
    
    size_t assignEnd = std::min(oldControl.numTrees, newControl.numTrees);
    
    // copy in trees that will be in the new one
    for (size_t treeNum = 0; treeNum < assignEnd; ++treeNum) {
      size_t oldTreeOffset = treeNum + oldSampleNum * oldControl.numTrees;
      size_t newTreeOffset = treeNum + newSampleNum * newControl.numTrees;
      
      // copy over children pointers, if any
      newTreeData.trees[newTreeOffset] = oldTreeData.trees[oldTreeOffset];
      setNewObservationIndices(newTreeData.trees[newTreeOffset].top,
                               newTreeData.treeIndices + newTreeOffset * data.numObservations,
                               oldTreeData.trees[oldTreeOffset].top);
      
      if (!newTreeData.trees[newTreeOffset].top.isBottom()) {
        newTreeData.trees[newTreeOffset].top.getRightChild()->parent = &newTreeData.trees[newTreeOffset].top;
        newTreeData.trees[newTreeOffset].top.getLeftChild()->parent  = &newTreeData.trees[newTreeOffset].top;
        
        // prevent destructor from freeing children since we just assigned them over
        oldTreeData.trees[oldTreeOffset].top.leftChild = NULL;
      }
    }
    
    size_t oldSampleOffset = oldSampleNum * data.numObservations * oldControl.numTrees;
    size_t newSampleOffset = newSampleNum * data.numObservations * newControl.numTrees;
    
    std::memcpy(newTreeData.treeIndices + newSampleOffset, oldTreeData.treeIndices + oldSampleOffset, assignEnd * data.numObservations * sizeof(size_t));
    std::memcpy(newTreeData.treeFits    + newSampleOffset, oldTreeData.treeFits    + oldSampleOffset, assignEnd * data.numObservations * sizeof(double));
    
    // if any new trees are required, create and initialize those
    for (size_t treeNum = assignEnd; treeNum < newControl.numTrees; ++treeNum) {
      size_t treeOffset = treeNum + newSampleNum * newControl.numTrees;
      new (newTreeData.trees + treeOffset)
        Tree(newTreeData.treeIndices + treeOffset * data.numObservations, data.numObservations, data.numPredictors);
      misc_setVectorToConstant(newTreeData.treeFits + treeOffset * data.numObservations, data.numObservations, 0.0);
    }
    
    // if any extra trees exist, delete them
    oldSampleOffset = oldSampleNum * oldControl.numTrees;
    
    for (size_t treeNum = oldControl.numTrees; treeNum > assignEnd; --treeNum)
      oldTreeData.trees[treeNum - 1 + oldSampleOffset].~Tree();
  }
}

namespace dbarts {
  bool State::resize(const BARTFit& fit, const Control& newControl) {
    const Control& oldControl(fit.control);
    const Data& data(fit.data);
   
    if (oldControl.keepTrees == newControl.keepTrees && oldControl.numTrees == newControl.numTrees) return false;
    
    State oldState = *this;
    
    if (oldControl.numTrees != newControl.numTrees) {
      treeIndices = new size_t[data.numObservations * newControl.numTrees];
      trees       = static_cast<Tree*>(::operator new (newControl.numTrees * sizeof(Tree)));
      treeFits    = new double[data.numObservations * newControl.numTrees];
      
      TreeData oldTrees = { oldState.treeIndices, oldState.trees, oldState.treeFits };
      TreeData newTrees = { treeIndices, trees, treeFits };
      ResizeData resizeData = { fit.data, oldControl, newControl, oldTrees, newTrees };
      
      copyTreesForSample(resizeData, 0, 0);
      
      delete [] oldState.treeFits;
      ::operator delete (oldState.trees);
      delete [] oldState.treeIndices;
    }
    
    if (newControl.keepTrees) {
      size_t totalNumTrees = newControl.numTrees * fit.currentNumSamples;
      
      savedTreeIndices = new size_t[data.numObservations * totalNumTrees];
      savedTrees = static_cast<Tree*>(::operator new (totalNumTrees * sizeof(Tree)));
      savedTreeFits = new double[data.numObservations * totalNumTrees];
      
      if (oldControl.keepTrees) {
        TreeData oldTrees = { oldState.savedTreeIndices, oldState.savedTrees, oldState.savedTreeFits };
        TreeData newTrees = { savedTreeIndices, savedTrees, savedTreeFits };
        ResizeData resizeData = { fit.data, oldControl, newControl, oldTrees, newTrees };
        
        for (size_t sampleNum = 0; sampleNum < fit.currentNumSamples; ++sampleNum)
          copyTreesForSample(resizeData, sampleNum, sampleNum);
        
        delete [] oldState.savedTreeFits;
        ::operator delete (oldState.savedTrees);
        delete [] oldState.savedTreeIndices;
      } else {
        for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum)
          new (savedTrees + treeNum) Tree(savedTreeIndices + treeNum * data.numObservations, data.numObservations, data.numPredictors);
        misc_setVectorToConstant(savedTreeFits, data.numObservations * totalNumTrees, 0.0);
      }
    } else {
      savedTreeIndices = NULL;
      savedTrees = NULL;
      savedTreeFits = NULL;
      
      if (oldControl.keepTrees) {
        delete [] oldState.savedTreeFits;
        for (size_t treeNum = oldControl.numTrees * fit.currentNumSamples; treeNum > 0; --treeNum)
          oldState.savedTrees[treeNum - 1].~Tree();
        ::operator delete (oldState.savedTrees);
        delete [] oldState.savedTreeIndices;
      }
    }
    
    return true;
  }
  
  bool State::resize(const BARTFit& fit, size_t newNumSamples) {
    const Control& control(fit.control);
    const Data& data(fit.data);
    
    size_t oldNumSamples = fit.currentNumSamples;
    
    if (newNumSamples == oldNumSamples || !control.keepTrees) return false;
    
    State oldState = *this;
    
    size_t totalNumTrees = control.numTrees * newNumSamples;
    
    savedTreeIndices = new size_t[data.numObservations * totalNumTrees];
    savedTrees = static_cast<Tree*>(::operator new (totalNumTrees * sizeof(Tree)));
    savedTreeFits = new double[data.numObservations * totalNumTrees];
    
    TreeData oldTrees = { oldState.savedTreeIndices, oldState.savedTrees, oldState.savedTreeFits };
    TreeData newTrees = { savedTreeIndices, savedTrees, savedTreeFits };
    ResizeData resizeData = { fit.data, control, control, oldTrees, newTrees };
        
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
    
    delete [] oldState.treeFits;
     ::operator delete (oldState.trees);
    delete [] oldState.treeIndices;
    
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
  const char* const* State::createTreeStrings(const BARTFit& fit, bool useSavedTrees) const
  {
    StringWriter writer;
    
    size_t numTrees;
    const Tree* targetTrees;
    if (!useSavedTrees || !fit.control.keepTrees) {
      numTrees = fit.control.numTrees;
      targetTrees = trees;
    } else {
      numTrees = fit.control.numTrees * fit.currentNumSamples;
      targetTrees = savedTrees;
    }
    
    char** result = new char*[numTrees];
    for (size_t i = 0; i < numTrees; ++i) {
      writer.buffer = new char[BASE_BUFFER_SIZE];
      writer.length = BASE_BUFFER_SIZE;
      writer.pos = 0;
      
      writer.writeNode(targetTrees[i].top);
     
      writer.writeChar('\0');
      
      result[i] = writer.buffer;
    }
    
    return(result);
  }
  
  void State::recreateTreesFromStrings(const BARTFit& fit, const char* const* treeStrings, bool useSavedTrees)
  {
    size_t numTrees;
    Tree* targetTrees;
    if (!useSavedTrees || !fit.control.keepTrees) {
      numTrees = fit.control.numTrees;
      targetTrees = trees;
    } else {
      numTrees = fit.control.numTrees * fit.currentNumSamples;
      targetTrees = savedTrees;
    }
    
    for (size_t i = 0; i < numTrees; ++i) {
      targetTrees[i].top.clear();
      readNode(targetTrees[i].top, treeStrings[i], fit.data.numPredictors);
      
      if (!targetTrees[i].top.isBottom()) {
        updateVariablesAvailable(fit, targetTrees[i].top, targetTrees[i].top.p.rule.variableIndex);
      
        targetTrees[i].top.addObservationsToChildren(fit);
      }
    }
  }
}
