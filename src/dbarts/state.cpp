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
#include <misc/memalign.h>
#include <misc/simd.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
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
    
    treeFitsAlignment = misc_simd_alignment;
    if (treeFitsAlignment == 0) {
      treeFitsStride = data.numObservations;
      treeFits = new double[treeFitsStride * totalNumTrees];
    } else {
      size_t remainder = data.numObservations % (treeFitsAlignment / sizeof(double));
      treeFitsStride = data.numObservations + 
        (remainder == 0 ? 0 : (treeFitsAlignment / sizeof(double) - remainder));
     if (misc_alignedAllocate(reinterpret_cast<void**>(&treeFits), treeFitsAlignment, treeFitsStride * totalNumTrees * sizeof(double)) != 0)
        ext_throwError("error allocating aligned vector");
    }
    misc_setVectorToConstant(treeFits, treeFitsStride * totalNumTrees, 0.0);
    
    if (control.keepTrees) {
      totalNumTrees *= control.defaultNumSamples;
      
      savedTrees = static_cast<SavedTree*>(::operator new (totalNumTrees * sizeof(SavedTree)));
      for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum)
        new (savedTrees + treeNum) SavedTree();
    } else {
      savedTrees = NULL;
    }
        
    rng = NULL;
  }
  
  void State::invalidate(size_t numTrees, size_t numSamples) {
    if (savedTrees != NULL) {
      for (size_t treeNum = numTrees * numSamples; treeNum > 0; --treeNum)
        savedTrees[treeNum - 1].~SavedTree();
      ::operator delete (savedTrees);
    }
    
    if (treeFitsAlignment == 0) {
      delete [] treeFits;
    } else {
      misc_alignedFree(treeFits);
    }
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
    
    size_t treeFitsStride;
  };
  struct SavedResizeData {
    const Data& data;
    
    const Control& oldControl;
    const Control& newControl;
    
    const SavedTree* oldTrees;
    SavedTree* newTrees;
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
    
    oldSampleOffset = oldSampleNum * resizeData.treeFitsStride * oldControl.numTrees;
    newSampleOffset = newSampleNum * resizeData.treeFitsStride * newControl.numTrees;
    std::memcpy(newTreeData.treeFits + newSampleOffset, oldTreeData.treeFits + oldSampleOffset, assignEnd * resizeData.treeFitsStride * sizeof(double));
    
    // if any new trees are required, create and initialize those
    for (size_t treeNum = assignEnd; treeNum < newControl.numTrees; ++treeNum) {
      size_t treeOffset = treeNum + newSampleNum * newControl.numTrees;
      new (newTreeData.trees + treeOffset)
        Tree(newTreeData.treeIndices + treeOffset * data.numObservations, data.numObservations, data.numPredictors);
      misc_setVectorToConstant(newTreeData.treeFits + treeOffset * resizeData.treeFitsStride, resizeData.treeFitsStride, 0.0);
    }
    
    // if any extra trees exist, delete them
    oldSampleOffset = oldSampleNum * oldControl.numTrees;
    
    for (size_t treeNum = oldControl.numTrees; treeNum > assignEnd; --treeNum)
      oldTreeData.trees[treeNum - 1 + oldSampleOffset].~Tree();
  }
  
  void copyTreesForSample(SavedResizeData& resizeData, size_t oldSampleNum, size_t newSampleNum) {
    const Control& oldControl(resizeData.oldControl);
    const Control& newControl(resizeData.newControl);
    
    const SavedTree* oldTrees = resizeData.oldTrees;
    SavedTree* newTrees       = resizeData.newTrees;
    
    size_t assignEnd = std::min(oldControl.numTrees, newControl.numTrees);
    
    // copy in trees that will be in the new one
    for (size_t treeNum = 0; treeNum < assignEnd; ++treeNum) {
      size_t oldTreeOffset = treeNum + oldSampleNum * oldControl.numTrees;
      size_t newTreeOffset = treeNum + newSampleNum * newControl.numTrees;
      
      // copy over children pointers, if any
      newTrees[newTreeOffset] = oldTrees[oldTreeOffset];
      
      if (!newTrees[newTreeOffset].top.isBottom()) {
        newTrees[newTreeOffset].top.getRightChild()->parent = &newTrees[newTreeOffset].top;
        newTrees[newTreeOffset].top.getLeftChild()->parent  = &newTrees[newTreeOffset].top;
        
        // prevent destructor from freeing children since we just assigned them over
        const_cast<SavedTree*>(oldTrees)[oldTreeOffset].top.leftChild = NULL;
      }
    }
    
    // if any new trees are required, create and initialize those
    for (size_t treeNum = assignEnd; treeNum < newControl.numTrees; ++treeNum) {
      size_t treeOffset = treeNum + newSampleNum * newControl.numTrees;
      new (newTrees + treeOffset) SavedTree();
    }
    
    // if any extra trees exist, delete them
    size_t oldSampleOffset = oldSampleNum * oldControl.numTrees;
    
    for (size_t treeNum = oldControl.numTrees; treeNum > assignEnd; --treeNum)
      oldTrees[treeNum - 1 + oldSampleOffset].~SavedTree();
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
      if (treeFitsAlignment == 0) {
        treeFits = new double[treeFitsStride * newControl.numTrees];
      } else {
        if (misc_alignedAllocate(
              reinterpret_cast<void**>(&treeFits),
              treeFitsAlignment,
              treeFitsStride * newControl.numTrees * sizeof(double)) != 0)
          ext_throwError("error allocating aligned vector");
      }
      TreeData oldTrees = { oldState.treeIndices, oldState.trees, oldState.treeFits };
      TreeData newTrees = { treeIndices, trees, treeFits };
      ResizeData resizeData = { fit.data, oldControl, newControl, oldTrees, newTrees, treeFitsStride };
      
      copyTreesForSample(resizeData, 0, 0);
      
      if (treeFitsAlignment == 0) {
        delete [] oldState.treeFits;
      } else {
        misc_alignedFree(oldState.treeFits);
      }
      ::operator delete (oldState.trees);
      delete [] oldState.treeIndices;
    }
    
    if (newControl.keepTrees) {
      size_t totalNumTrees = newControl.numTrees * fit.currentNumSamples;
      
      savedTrees = static_cast<SavedTree*>(::operator new (totalNumTrees * sizeof(SavedTree)));
      
      if (oldControl.keepTrees) {
        SavedResizeData resizeData = { fit.data, oldControl, newControl, oldState.savedTrees, savedTrees };
        
        for (size_t sampleNum = 0; sampleNum < fit.currentNumSamples; ++sampleNum)
          copyTreesForSample(resizeData, sampleNum, sampleNum);
        
        ::operator delete (oldState.savedTrees);
      } else {
        for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum)
          new (savedTrees + treeNum) SavedTree();
      }
    } else {
      savedTrees = NULL;
      
      if (oldControl.keepTrees) {
        for (size_t treeNum = oldControl.numTrees * fit.currentNumSamples; treeNum > 0; --treeNum)
          oldState.savedTrees[treeNum - 1].~SavedTree();
        ::operator delete (oldState.savedTrees);
      }
    }
    
    return true;
  }
  
  bool State::resize(const BARTFit& fit, size_t newNumSamples) {
    const Control& control(fit.control);
    
    size_t oldNumSamples = fit.currentNumSamples;
    
    if (newNumSamples == oldNumSamples || !control.keepTrees) return false;
    
    State oldState = *this;
    
    size_t totalNumTrees = control.numTrees * newNumSamples;
    
    savedTrees = static_cast<SavedTree*>(::operator new (totalNumTrees * sizeof(SavedTree)));
    
    SavedResizeData resizeData = { fit.data, control, control, oldState.savedTrees, savedTrees };
        
    size_t numSamplesToCopy, oldSampleStart, newSampleStart;
    if (oldNumSamples > newNumSamples) {
      numSamplesToCopy = newNumSamples;
      oldSampleStart = oldNumSamples - newNumSamples;
      newSampleStart = 0;
      
      for (size_t sampleNum = oldNumSamples - newNumSamples; sampleNum > 0; --sampleNum) {
        for (size_t treeNum = control.numTrees; treeNum > 0; --treeNum) {
          size_t treeOffset = (treeNum - 1) + (sampleNum - 1) * control.numTrees;
          oldState.savedTrees[treeOffset].~SavedTree();
        }
      }
    } else {
      numSamplesToCopy = oldNumSamples;
      oldSampleStart = 0;
      newSampleStart = newNumSamples - oldNumSamples;
      
      for (size_t sampleNum = 0; sampleNum < newNumSamples - oldNumSamples; ++sampleNum) {
        size_t sampleOffset = sampleNum * control.numTrees;
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          new (savedTrees + treeNum + sampleOffset) SavedTree();
        }
      }
    }
        
    for (size_t sampleNum = 0; sampleNum < numSamplesToCopy; ++sampleNum)
      copyTreesForSample(resizeData, oldSampleStart + sampleNum, newSampleStart + sampleNum);
    
    ::operator delete (oldState.savedTrees);
    
    return true;
  }
  
  size_t State::getSerializedTreesLength(const BARTFit& fit) const {
    const Control& control(fit.control);
        
    size_t stateLength = 0;
    
    for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum)
      stateLength += trees[treeNum].getSerializedLength(fit);
    
    return stateLength;
  }
  
  void State::serializeTrees(const BARTFit& fit, void* state) const {
    const Control& control(fit.control);
    
    for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
      size_t stateLength = trees[treeNum].serialize(fit, state);
      state = reinterpret_cast<void*>(reinterpret_cast<char*>(state) + stateLength);
    }
  }
  
  void State::deserializeTrees(const BARTFit& fit, const void* state) {
    const Control& control(fit.control);
    
    for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
      size_t stateLength = trees[treeNum].deserialize(fit, state);
      state = reinterpret_cast<const void*>(reinterpret_cast<const char*>(state) + stateLength);
    }
  }
  
  size_t State::getSerializedSavedTreesLength(const BARTFit& fit) const {
    const Control& control(fit.control);
        
    size_t stateLength = 0;
           
    if (control.keepTrees) {
      size_t numSavedTrees = fit.currentNumSamples * control.numTrees;
      for (size_t treeNum = 0; treeNum < numSavedTrees; ++treeNum) {
        stateLength += savedTrees[treeNum].getSerializedLength();
      }
    }
    
    return stateLength;
  }
  
  void State::serializeSavedTrees(const BARTFit& fit, void* state) const {
    const Control& control(fit.control);
    
    if (control.keepTrees) {
      size_t numSavedTrees = fit.currentNumSamples * control.numTrees;
      for (size_t treeNum = 0; treeNum < numSavedTrees; ++treeNum) {
        size_t stateLength = savedTrees[treeNum].serialize(state);
        state = reinterpret_cast<void*>(reinterpret_cast<char*>(state) + stateLength);
      }
    }
  }
  
  void State::deserializeSavedTrees(const BARTFit& fit, const void* state) {
    const Control& control(fit.control);
    
    if (control.keepTrees) {
      size_t numSavedTrees = fit.currentNumSamples * control.numTrees;
      for (size_t treeNum = 0; treeNum < numSavedTrees; ++treeNum) {
        size_t stateLength = savedTrees[treeNum].deserialize(state);
        state = reinterpret_cast<const void*>(reinterpret_cast<const char*>(state) + stateLength);
      }
    }
  }
}
