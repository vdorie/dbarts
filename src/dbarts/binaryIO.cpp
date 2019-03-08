#include "binaryIO.hpp"

#include <cstddef>
#include <dbarts/cstdint.hpp>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include <misc/alloca.h>
#include <misc/binaryIO.h>

#include <external/io.h>
#include <external/random.h>

#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/state.hpp>

#include "tree.hpp"

#define VERSION_STRING_LENGTH 8

using std::size_t;

namespace dbarts {
#define CONTROL_BINARY_RESPONSE 1
#define CONTROL_VERBOSE         2
#define CONTROL_KEEP_TRAINING   4
#define CONTROL_USE_QUANTILES   8
  
  bool writeControl(misc_binaryIO* bio, const Control& control) {
    int errorCode = 0;
    
    uint32_t controlFlags = 0;
    controlFlags += control.responseIsBinary ? CONTROL_BINARY_RESPONSE : 0;
    controlFlags += control.verbose ? CONTROL_VERBOSE : 0;
    controlFlags += control.keepTrainingFits ? CONTROL_KEEP_TRAINING : 0;
    controlFlags += control.useQuantiles ? CONTROL_USE_QUANTILES : 0;
    
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, controlFlags)) != 0) goto write_control_cleanup;
    
    if ((errorCode = misc_bio_writeSizeType(bio, control.defaultNumSamples)) != 0) goto write_control_cleanup;
    if ((errorCode = misc_bio_writeSizeType(bio, control.defaultNumBurnIn)) != 0) goto write_control_cleanup;
    if ((errorCode = misc_bio_writeSizeType(bio, control.numTrees)) != 0) goto write_control_cleanup;
    
    // new as of version 00.09.00
    if ((errorCode = misc_bio_writeSizeType(bio, control.numChains)) != 0) goto write_control_cleanup;
    
    
    if ((errorCode = misc_bio_writeSizeType(bio, control.numThreads)) != 0) goto write_control_cleanup;
    
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, control.treeThinningRate)) != 0) goto write_control_cleanup;
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, control.printEvery)) != 0) goto write_control_cleanup;
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, control.printCutoffs)) != 0) goto write_control_cleanup;
    
    // new as of version 00.09.00
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, static_cast<uint32_t>(control.rng_algorithm))) != 0) goto write_control_cleanup;
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, static_cast<uint32_t>(control.rng_standardNormal))) != 0) goto write_control_cleanup;
    
    
write_control_cleanup:
    
    if (errorCode != 0) ext_issueWarning("error writing control object: %s", std::strerror(errorCode));
    
    return errorCode == 0;
  }
  
  bool readControl(misc_binaryIO* bio, Control& control, const Version& version) {
    int errorCode = 0;
    uint32_t enumTemp;
    
    uint32_t controlFlags = 0;
    controlFlags += control.responseIsBinary ? CONTROL_BINARY_RESPONSE : 0;
    controlFlags += control.verbose ? CONTROL_VERBOSE : 0;
    controlFlags += control.keepTrainingFits ? CONTROL_KEEP_TRAINING : 0;
    controlFlags += control.useQuantiles ? CONTROL_USE_QUANTILES : 0;
    
    if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, &controlFlags)) != 0) goto read_control_cleanup;
    control.responseIsBinary = (controlFlags & CONTROL_BINARY_RESPONSE) != 0;
    control.verbose = (controlFlags & CONTROL_VERBOSE) != 0;
    control.keepTrainingFits = (controlFlags & CONTROL_KEEP_TRAINING) != 0;
    control.useQuantiles = (controlFlags & CONTROL_USE_QUANTILES) != 0;
    
    if ((errorCode = misc_bio_readSizeType(bio, &control.defaultNumSamples)) != 0) goto read_control_cleanup;
    if ((errorCode = misc_bio_readSizeType(bio, &control.defaultNumBurnIn)) != 0) goto read_control_cleanup;
    if ((errorCode = misc_bio_readSizeType(bio, &control.numTrees)) != 0) goto read_control_cleanup;
    
    if (version.major > 0 || version.minor > 8) {
      if ((errorCode = misc_bio_readSizeType(bio, &control.numChains)) != 0) goto read_control_cleanup;
    } else {
      control.numChains = 1;
    }
    
    if ((errorCode = misc_bio_readSizeType(bio, &control.numThreads)) != 0) goto read_control_cleanup;
    
    if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, &control.treeThinningRate)) != 0) goto read_control_cleanup;
    if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, &control.printEvery)) != 0) goto read_control_cleanup;
    if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, &control.printCutoffs)) != 0) goto read_control_cleanup;
    
    if (version.major > 0 || version.minor > 8) {
      if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, &enumTemp)) != 0) goto read_control_cleanup;
      control.rng_algorithm = static_cast<rng_algorithm_t>(enumTemp);
      if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, &enumTemp)) != 0) goto read_control_cleanup;
      control.rng_standardNormal = static_cast<rng_standardNormal_t>(enumTemp);
    } else {
      control.rng_algorithm = RNG_ALGORITHM_INVALID;
      control.rng_standardNormal = RNG_STANDARD_NORMAL_INVALID;
    }
    
    control.callback = NULL;
    control.callbackData = NULL;
    
read_control_cleanup:
    
    if (errorCode != 0) ext_issueWarning("error reading control object: %s", std::strerror(errorCode));

    return errorCode == 0;
  }
  
#define DATA_HAS_WEIGHTS      1
#define DATA_HAS_OFFSET       2
#define DATA_HAS_TEST_OFFSET  4
#define DATA_HAS_MAX_NUM_CUTS 8
    
  bool writeData(misc_binaryIO* bio, const Data& data) {
    int errorCode = 0;
    uint32_t* variableTypes = NULL;
    
    uint32_t dataFlags = 0;
    dataFlags |= ((data.weights != NULL) ? DATA_HAS_WEIGHTS : 0);
    dataFlags |= ((data.offset != NULL) ?  DATA_HAS_OFFSET : 0);
    dataFlags |= ((data.testOffset != NULL) ? DATA_HAS_TEST_OFFSET : 0);
    dataFlags |= ((data.maxNumCuts != NULL) ? DATA_HAS_MAX_NUM_CUTS : 0);
    
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, dataFlags)) != 0) goto write_data_cleanup;
    
    if ((errorCode = misc_bio_writeSizeType(bio, data.numObservations)) != 0) goto write_data_cleanup;
    if ((errorCode = misc_bio_writeSizeType(bio, data.numPredictors)) != 0) goto write_data_cleanup;
    if ((errorCode = misc_bio_writeSizeType(bio, data.numTestObservations)) != 0) goto write_data_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, data.sigmaEstimate)) != 0) goto write_data_cleanup;
    
    if ((errorCode = misc_bio_writeNDoubles(bio, data.y, data.numObservations)) != 0) goto write_data_cleanup;
    if ((errorCode = misc_bio_writeNDoubles(bio, data.x, data.numObservations * data.numPredictors)) != 0) goto write_data_cleanup;
    if (data.numTestObservations > 0 &&
      (errorCode = misc_bio_writeNDoubles(bio, data.x_test, data.numTestObservations * data.numPredictors)) != 0) goto write_data_cleanup;
    
    if (data.weights != NULL && 
      (errorCode = misc_bio_writeNDoubles(bio, data.weights, data.numObservations)) != 0) goto write_data_cleanup;
    if (data.offset != NULL && 
      (errorCode = misc_bio_writeNDoubles(bio, data.offset, data.numObservations)) != 0) goto write_data_cleanup;
    if (data.testOffset != NULL && 
      (errorCode = misc_bio_writeNDoubles(bio, data.testOffset, data.numTestObservations)) != 0) goto write_data_cleanup;
    
    variableTypes = misc_stackAllocate(data.numPredictors, uint32_t);
    for (size_t j = 0; j < data.numPredictors; ++j) variableTypes[j] = static_cast<uint32_t>(data.variableTypes[j]);
    
    if ((errorCode = misc_bio_writeNUnsigned32BitIntegers(bio, variableTypes, data.numPredictors)) != 0) goto write_data_cleanup;
    
    errorCode = misc_bio_writeNUnsigned32BitIntegers(bio, data.maxNumCuts, data.numPredictors);
    
write_data_cleanup:
    if (variableTypes != NULL) { misc_stackFree(variableTypes); }

    if (errorCode != 0) ext_issueWarning("error writing data object: %s", std::strerror(errorCode));
    
    return errorCode == 0;
  }
  
  bool readData(misc_binaryIO* bio, Data& data) {
    int errorCode = 0;
    uint32_t* variableTypes = NULL;
    
    uint32_t dataFlags = 0;
    
    if ((errorCode = misc_bio_readNUnsigned32BitIntegers(bio, &dataFlags, 1)) != 0) goto read_data_cleanup;
    
    if ((errorCode = misc_bio_readSizeType(bio, &data.numObservations)) != 0) goto read_data_cleanup;
    if ((errorCode = misc_bio_readSizeType(bio, &data.numPredictors)) != 0) goto read_data_cleanup;
    if ((errorCode = misc_bio_readSizeType(bio, &data.numTestObservations)) != 0) goto read_data_cleanup;
    if ((errorCode = misc_bio_readDouble(bio, &data.sigmaEstimate)) != 0) goto read_data_cleanup;
    
    
    data.y = new double[data.numObservations];
    if ((errorCode = misc_bio_readNDoubles(bio, const_cast<double*>(data.y), data.numObservations)) != 0) goto read_data_cleanup;
    
    data.x = new double[data.numObservations * data.numPredictors];
    if ((errorCode = misc_bio_readNDoubles(bio, const_cast<double*>(data.x), data.numObservations * data.numPredictors)) != 0) goto read_data_cleanup;
    
    if (data.numTestObservations > 0) {
      data.x_test = new double[data.numTestObservations * data.numPredictors];
      if ((errorCode = misc_bio_readNDoubles(bio, const_cast<double*>(data.x_test), data.numTestObservations * data.numPredictors)) != 0) goto read_data_cleanup;
    } else data.x_test = NULL;
    
    if (dataFlags & DATA_HAS_WEIGHTS) {
      data.weights = new double[data.numObservations];
      if ((errorCode = misc_bio_readNDoubles(bio, const_cast<double*>(data.weights), data.numObservations)) != 0) goto read_data_cleanup;
    } else data.weights = NULL;
    
    if (dataFlags & DATA_HAS_OFFSET) {
      data.offset = new double[data.numObservations];
      if ((errorCode = misc_bio_readNDoubles(bio, const_cast<double*>(data.offset), data.numObservations)) != 0) goto read_data_cleanup;
    } else data.offset = NULL;
    
    if (dataFlags & DATA_HAS_TEST_OFFSET) {
      data.testOffset = new double[data.numTestObservations];
      if ((errorCode = misc_bio_readNDoubles(bio, const_cast<double*>(data.testOffset), data.numTestObservations)) != 0) goto read_data_cleanup;
    } else data.testOffset = NULL;
    
    variableTypes = misc_stackAllocate(data.numPredictors, uint32_t);
    if ((errorCode = misc_bio_readNUnsigned32BitIntegers(bio, variableTypes, data.numPredictors)) != 0) goto read_data_cleanup;
    data.variableTypes = new VariableType[data.numPredictors];
    for (size_t j = 0; j < data.numPredictors; ++j) const_cast<VariableType*>(data.variableTypes)[j] = static_cast<VariableType>(variableTypes[j]);
    
    if (dataFlags & DATA_HAS_MAX_NUM_CUTS) {
      data.maxNumCuts = new uint32_t[data.numPredictors];
      if ((errorCode = misc_bio_readNUnsigned32BitIntegers(bio, const_cast<uint32_t*>(data.maxNumCuts), data.numPredictors)) != 0) goto read_data_cleanup; 
    } else data.maxNumCuts = NULL;
    
read_data_cleanup:
    if (variableTypes != NULL) { misc_stackFree(variableTypes); }
      
    if (errorCode != 0) {
      delete [] data.maxNumCuts;
      delete [] data.variableTypes;
      delete [] data.testOffset;
      delete [] data.offset;
      delete [] data.weights;
      delete [] data.x_test;
      delete [] data.x;
      delete [] data.y;
    
      ext_issueWarning("error reading data object: %s", std::strerror(errorCode));
    }
    
    return errorCode == 0;
  }
  
  bool writeModel(misc_binaryIO* bio, const Model& model)
  {
    int errorCode = 0;
    
    if ((errorCode = misc_bio_writeDouble(bio, model.birthOrDeathProbability)) != 0) goto write_model_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, model.swapProbability)) != 0) goto write_model_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, model.changeProbability)) != 0) goto write_model_cleanup;
    
    if ((errorCode = misc_bio_writeDouble(bio, model.birthProbability)) != 0) goto write_model_cleanup;
    
    // this needs some seeerious work
    if ((errorCode = misc_bio_writeNChars(bio, "cgm ", 4)) != 0) goto write_model_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, static_cast<CGMPrior*>(model.treePrior)->base)) != 0) goto write_model_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, static_cast<CGMPrior*>(model.treePrior)->power)) != 0) goto write_model_cleanup;
    
    
    if ((errorCode = misc_bio_writeNChars(bio, "nrml", 4)) != 0) goto write_model_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, static_cast<NormalPrior*>(model.muPrior)->precision)) != 0) goto write_model_cleanup;
    
    if ((errorCode = misc_bio_writeNChars(bio, "chsq", 4)) != 0) goto write_model_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, static_cast<ChiSquaredPrior*>(model.sigmaSqPrior)->degreesOfFreedom)) != 0) goto write_model_cleanup;
    if ((errorCode = misc_bio_writeDouble(bio, static_cast<ChiSquaredPrior*>(model.sigmaSqPrior)->scale)) != 0) goto write_model_cleanup;
    
write_model_cleanup:
    if (errorCode != 0) ext_issueWarning("error writing model object: %s", std::strerror(errorCode));
    
    return errorCode == 0;
  }
  
  bool readModel(misc_binaryIO* bio, Model& model)
  {
    int errorCode = 0;
    char priorName[4];
    
    if ((errorCode = misc_bio_readDouble(bio, &model.birthOrDeathProbability)) != 0) goto read_model_cleanup;
    if ((errorCode = misc_bio_readDouble(bio, &model.swapProbability)) != 0) goto read_model_cleanup;
    if ((errorCode = misc_bio_readDouble(bio, &model.changeProbability)) != 0) goto read_model_cleanup;
    
    if ((errorCode = misc_bio_readDouble(bio, &model.birthProbability)) != 0) goto read_model_cleanup;
    
    if ((errorCode = misc_bio_readNChars(bio, priorName, 4)) != 0) goto read_model_cleanup;
    if (std::strncmp(priorName, "cgm ", 4) != 0) { errorCode = EILSEQ; goto read_model_cleanup; }
    
    model.treePrior = new CGMPrior;
    if ((errorCode = misc_bio_readDouble(bio, &static_cast<CGMPrior*>(model.treePrior)->base)) != 0) goto read_model_cleanup;
    if ((errorCode = misc_bio_readDouble(bio, &static_cast<CGMPrior*>(model.treePrior)->power)) != 0) goto read_model_cleanup;
    
    
    if ((errorCode = misc_bio_readNChars(bio, priorName, 4)) != 0) goto read_model_cleanup;
    if (std::strncmp(priorName, "nrml", 4) != 0) { errorCode = EILSEQ; goto read_model_cleanup; }
    
    model.muPrior = new NormalPrior;
    if ((errorCode = misc_bio_readDouble(bio, &static_cast<NormalPrior*>(model.muPrior)->precision)) != 0) goto read_model_cleanup;
    
    
    if ((errorCode = misc_bio_readNChars(bio, priorName, 4)) != 0) goto read_model_cleanup;
    if (std::strncmp(priorName, "chsq", 4) != 0) { errorCode = EILSEQ; goto read_model_cleanup; }
    
    model.sigmaSqPrior = new ChiSquaredPrior;
    if ((errorCode = misc_bio_readDouble(bio, &static_cast<ChiSquaredPrior*>(model.sigmaSqPrior)->degreesOfFreedom)) != 0) goto read_model_cleanup;
    if ((errorCode = misc_bio_readDouble(bio, &static_cast<ChiSquaredPrior*>(model.sigmaSqPrior)->scale)) != 0) goto read_model_cleanup;
    
read_model_cleanup:
    
    if (errorCode != 0) {
      delete model.sigmaSqPrior;
      delete model.muPrior;
      delete model.treePrior;
      
      ext_issueWarning("error reading model object: %s", std::strerror(errorCode));
    }
    
    return errorCode == 0;
  }
}

namespace {
  int writeTree(misc_binaryIO* bio, const dbarts::Tree& tree, const dbarts::Data& data, const size_t* treeIndices);
  int readTree(misc_binaryIO* bio, dbarts::Tree& tree, const dbarts::Data& data, const size_t* treeIndices);
  int writeTree(misc_binaryIO* bio, const dbarts::SavedTree& tree);
  int readTree(misc_binaryIO* bio, dbarts::SavedTree& tree);
}

namespace dbarts {
  bool writeState(misc_binaryIO* bio, const State* state, const Control& control, const Data& data, size_t numSamples)
  {
    int errorCode = 0;
    size_t treeNum;
    size_t rngStateLength;
    int* rngState = NULL;
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      if ((errorCode = misc_bio_writeNSizeTypes(bio, state[chainNum].treeIndices, data.numObservations * control.numTrees)) != 0) goto write_state_cleanup;
      for (treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        if ((errorCode = writeTree(bio, state[chainNum].trees[treeNum], data, state[chainNum].treeIndices + treeNum * data.numObservations)) != 0) goto write_state_cleanup;
      }
      if ((errorCode = misc_bio_writeNDoubles(bio, state[chainNum].treeFits, data.numObservations * control.numTrees)) != 0) goto write_state_cleanup;
      
      // saved trees
      if (control.keepTrees) {
        size_t totalNumTrees = numSamples * control.numTrees;
        for (treeNum = 0; treeNum < totalNumTrees; ++treeNum)
          if ((errorCode = writeTree(bio, state[chainNum].savedTrees[treeNum])) != 0) goto write_state_cleanup;
      }
      
      
      if ((errorCode = misc_bio_writeDouble(bio, state[chainNum].sigma)) != 0) goto write_state_cleanup;
      
      rngStateLength = ext_rng_getSerializedStateLength(state[chainNum].rng) / sizeof(int);
      if ((errorCode = misc_bio_writeSizeType(bio, rngStateLength)) != 0) goto write_state_cleanup;
      
      rngState = new int[rngStateLength];
      if (rngState == NULL) {
        errorCode = ENOMEM;
        goto write_state_cleanup;
      }
      ext_rng_writeSerializedState(state[chainNum].rng, static_cast<void*>(rngState));
      
      if ((errorCode = misc_bio_writeNInts(bio, rngState, rngStateLength)) != 0) {
        delete [] rngState;
        goto write_state_cleanup;
      }
      
      delete [] rngState;
      rngState = NULL;
    }
    
write_state_cleanup:
    if (errorCode != 0) ext_issueWarning("error writing state object: %s", std::strerror(errorCode));
    
    return errorCode == 0;
  }
  
  bool readState(misc_binaryIO* bio, State* state, const Control& control, const Data& data, size_t numSamples, const Version& version)
  {
    int errorCode = 0;
    size_t treeNum;
    size_t rngStateLength;
    int* rngState = NULL;
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      if ((errorCode = misc_bio_readNSizeTypes(bio, state[chainNum].treeIndices, data.numObservations * control.numTrees)) != 0) goto read_state_cleanup;
      
      for (treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        if ((errorCode = readTree(bio, state[chainNum].trees[treeNum], data, state[chainNum].treeIndices + treeNum * data.numObservations)) != 0) goto read_state_cleanup;
      }
      if ((errorCode = misc_bio_readNDoubles(bio, state[chainNum].treeFits, data.numObservations * control.numTrees)) != 0) goto read_state_cleanup;
      
      if (control.keepTrees) {
        size_t totalNumTrees = numSamples * control.numTrees;
        for (treeNum = 0; treeNum < totalNumTrees; ++treeNum)
          if ((errorCode = readTree(bio, state[chainNum].savedTrees[treeNum])) != 0) goto read_state_cleanup;
      }
      
      if ((errorCode = misc_bio_readDouble(bio, &state[chainNum].sigma)) != 0) goto read_state_cleanup;
      
      if (version.major > 0 || version.minor > 8) {
        if ((errorCode = misc_bio_readSizeType(bio, &rngStateLength)) != 0) goto read_state_cleanup;
        
        rngState = new int[rngStateLength];
        if (rngState == NULL) {
          errorCode = ENOMEM;
          goto read_state_cleanup;
        }
        if ((errorCode = misc_bio_readNInts(bio, rngState, rngStateLength)) != 0) {
          delete [] rngState;
          goto read_state_cleanup;
        }
        
        ext_rng_readSerializedState(state[chainNum].rng, static_cast<const void*>(rngState));
        
        delete [] rngState;
        rngState = NULL;
      }
    }

read_state_cleanup:
    if (errorCode != 0) ext_issueWarning("error reading state object: %s", std::strerror(errorCode));
    
    return errorCode == 0;
  }
  
  // version should be of format MM.mm.rr
  //   MM - major
  //   mm - minor
  //   rr - revision
  int readVersion(misc_binaryIO* bio, Version& version, char** versionStringPtr)
  {
    char versionString[VERSION_STRING_LENGTH + 1];
    versionString[VERSION_STRING_LENGTH] = '\0';
    
    char* originalVersionString = new char[VERSION_STRING_LENGTH + 1];
    for (size_t i = 0; i < VERSION_STRING_LENGTH + 1; ++i)
     originalVersionString[i] = '\0'; 
    
    int errorCode;
    if ((errorCode = misc_bio_readNChars(bio, versionString, VERSION_STRING_LENGTH)) != 0)
      goto misc_bio_readVersion_error;
    
    std::memcpy(originalVersionString, const_cast<const char*>(versionString), (VERSION_STRING_LENGTH + 1) * sizeof(char));
    
    
    if (versionString[2] != '.' || versionString[5] != '.') {
      errorCode = EINVAL;
      goto misc_bio_readVersion_error;
    }
    
    versionString[2] = '\0';
    versionString[5] = '\0';
    
    errno = 0;
    version.major = static_cast<size_t>(std::strtol(versionString, NULL, 10));
    if (version.major == 0 && errno != 0) {
      errorCode = errno;
      goto misc_bio_readVersion_error;
    }
    
    version.minor = static_cast<size_t>(std::strtol(versionString + 3, NULL, 10));
    if (version.minor == 0 && errno != 0) {
      errorCode = errno;
      goto misc_bio_readVersion_error;
    }
    
    version.revision = static_cast<size_t>(std::strtol(versionString + 6, NULL, 10));
    if (version.revision == 0 && errno != 0) {
      errorCode = errno;
      goto misc_bio_readVersion_error;
    }
    
    if (version.major == 0 && version.minor == 0 && version.revision == 0) {
      errorCode = EDOM;
      goto misc_bio_readVersion_error;
    }
    
    delete [] originalVersionString;
    
    return 0;

misc_bio_readVersion_error:
    if (versionStringPtr != NULL) {
      *versionStringPtr = originalVersionString;
    } else {
      delete [] originalVersionString;
    }
    
    return errorCode;
  }
}

#define NODE_HAS_CHILDREN 1
namespace {
  int writeNode(misc_binaryIO* bio, const dbarts::Node& node, const dbarts::Data& data, const size_t* treeIndices)
  {
    int errorCode = 0;
    
    ptrdiff_t observationOffset = 0;
    unsigned char nodeFlags = 0;
    uint64_t variablesAvailableForSplit = 0;
    
    observationOffset = node.observationIndices - treeIndices;
    if (observationOffset < 0) {
      errorCode = EINVAL; goto write_node_cleanup;
    } else if ((errorCode = misc_bio_writeSizeType(bio, static_cast<size_t>(observationOffset))) != 0) goto write_node_cleanup;
    
    if ((errorCode = misc_bio_writeSizeType(bio, node.enumerationIndex)) != 0) goto write_node_cleanup;
    if ((errorCode = misc_bio_writeSizeType(bio, node.numObservations)) != 0) goto write_node_cleanup;
    
    for (size_t j = 0; j < data.numPredictors; ++j) {
      if (node.variablesAvailableForSplit[j] == true) variablesAvailableForSplit |= 1 << j;
    }
    if ((errorCode = misc_bio_writeUnsigned64BitInteger(bio, variablesAvailableForSplit)) != 0) goto write_node_cleanup;
    
    if (node.leftChild != NULL) {
      nodeFlags += NODE_HAS_CHILDREN;
      
      if ((errorCode = misc_bio_writeChar(bio, *reinterpret_cast<char*>(&nodeFlags))) != 0) goto write_node_cleanup;
      
      if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, *(reinterpret_cast<const uint32_t*>(&node.p.rule.variableIndex)))) != 0) goto write_node_cleanup;
      if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, node.p.rule.categoryDirections)) != 0) goto write_node_cleanup;
      
      if ((errorCode = writeNode(bio, *node.leftChild, data, treeIndices))) goto write_node_cleanup;
      if ((errorCode = writeNode(bio, *node.p.rightChild, data, treeIndices))) goto write_node_cleanup;
    } else {
      if ((errorCode = misc_bio_writeChar(bio, *reinterpret_cast<char*>(&nodeFlags))) != 0) goto write_node_cleanup;
      
      if ((errorCode = misc_bio_writeDouble(bio, node.m.average)) != 0) goto write_node_cleanup;
      if ((errorCode = misc_bio_writeDouble(bio, node.m.numEffectiveObservations)) != 0) goto write_node_cleanup;
    }
    
write_node_cleanup:
    
    return errorCode;
  }
  
  int writeNode(misc_binaryIO* bio, const dbarts::SavedNode& node)
  {
    int errorCode = 0;
    
    unsigned char nodeFlags = 0;
    
    if (node.getLeftChild() != NULL)
      nodeFlags += NODE_HAS_CHILDREN;
    if ((errorCode = misc_bio_writeChar(bio, *reinterpret_cast<char*>(&nodeFlags))) != 0) goto write_saved_node_cleanup;
    
    if ((errorCode = misc_bio_writeUnsigned32BitInteger(bio, *(reinterpret_cast<const uint32_t*>(&node.variableIndex)))) != 0) goto write_saved_node_cleanup;
    
    if ((errorCode = misc_bio_writeDouble(bio, node.split)) != 0) goto write_saved_node_cleanup;
    
    if (node.getLeftChild() != NULL) {
      if ((errorCode = writeNode(bio, *node.getLeftChild()))) goto write_saved_node_cleanup;
      if ((errorCode = writeNode(bio, *node.getRightChild()))) goto write_saved_node_cleanup;
    }
    
write_saved_node_cleanup:
    
    return errorCode;
  }
  
  int writeTree(misc_binaryIO* bio, const dbarts::Tree& tree, const dbarts::Data& data, const size_t* treeIndices)
  {
    return writeNode(bio, tree.top, data, treeIndices);
  }
  
  int writeTree(misc_binaryIO* bio, const dbarts::SavedTree& tree)
  {
    return writeNode(bio, tree.top);
  }
  
  int readNode(misc_binaryIO* bio, dbarts::Node& node, const dbarts::Data& data, const size_t* treeIndices)
  {
    int errorCode = 0;
    
    size_t observationOffset = 0;
    unsigned char nodeFlags = 0;
    uint64_t variablesAvailableForSplit = 0;
    dbarts::Node* leftChild = NULL;
    dbarts::Node* rightChild = NULL;
    
    if ((errorCode = misc_bio_readSizeType(bio, &observationOffset)) != 0) goto read_node_cleanup;
    if (observationOffset >= data.numObservations) { errorCode = EINVAL; goto read_node_cleanup; }
    node.observationIndices = const_cast<size_t*>(treeIndices) + observationOffset;
    
    if ((errorCode = misc_bio_readSizeType(bio, &node.enumerationIndex)) != 0) goto read_node_cleanup;
    if ((errorCode = misc_bio_readSizeType(bio, &node.numObservations)) != 0) goto read_node_cleanup;
    
    if ((errorCode = misc_bio_readUnsigned64BitInteger(bio, &variablesAvailableForSplit)) != 0) goto read_node_cleanup;
    for (size_t j = 0; j < data.numPredictors; ++j) {
      node.variablesAvailableForSplit[j] = (variablesAvailableForSplit & (1 << j)) != 0;
    }
    
    if ((errorCode = misc_bio_readChar(bio, reinterpret_cast<char*>(&nodeFlags))) != 0) goto read_node_cleanup;
    
    if (nodeFlags > NODE_HAS_CHILDREN) { errorCode = EINVAL; goto read_node_cleanup; }
    
    if (nodeFlags & NODE_HAS_CHILDREN) {
      if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, reinterpret_cast<uint32_t*>(&node.p.rule.variableIndex))) != 0) goto read_node_cleanup;
      if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, &node.p.rule.categoryDirections)) != 0) goto read_node_cleanup;
      
      leftChild = new dbarts::Node(node, data.numPredictors);
      node.leftChild = leftChild;
      if ((errorCode = readNode(bio, *leftChild, data, treeIndices)) != 0) goto read_node_cleanup;
      
      rightChild = new dbarts::Node(node, data.numPredictors);
      node.p.rightChild = rightChild;
      if ((errorCode = readNode(bio, *rightChild, data, treeIndices)) != 0) goto read_node_cleanup;
    } else {
      node.leftChild = NULL;
      if ((errorCode = misc_bio_readDouble(bio, &node.m.average)) != 0) goto read_node_cleanup;
      if ((errorCode = misc_bio_readDouble(bio, &node.m.numEffectiveObservations)) != 0) goto read_node_cleanup;  
    }
    
read_node_cleanup:

    if (errorCode != 0) {
      delete rightChild;
      delete leftChild;
      
      node.leftChild = NULL;
    }
    
    return errorCode;
  }
  
  int readNode(misc_binaryIO* bio, dbarts::SavedNode& node)
  {
    int errorCode = 0;
    
    unsigned char nodeFlags = 0;
    
    node.leftChild = NULL;
    node.rightChild = NULL;
    
    if ((errorCode = misc_bio_readChar(bio, reinterpret_cast<char*>(&nodeFlags))) != 0) goto read_saved_node_cleanup;
    if (nodeFlags > NODE_HAS_CHILDREN) { errorCode = EINVAL; goto read_saved_node_cleanup; }
    
    if ((errorCode = misc_bio_readUnsigned32BitInteger(bio, (reinterpret_cast<uint32_t*>(&node.variableIndex)))) != 0) goto read_saved_node_cleanup;
    
    if ((errorCode = misc_bio_readDouble(bio, &node.split)) != 0) goto read_saved_node_cleanup;
    
    if (nodeFlags & NODE_HAS_CHILDREN) {
      node.leftChild = new dbarts::SavedNode();
      node.leftChild->parent = &node;
      if ((errorCode = writeNode(bio, *node.getLeftChild()))) goto read_saved_node_cleanup;
      node.rightChild = new dbarts::SavedNode();
      node.rightChild->parent = &node;
      if ((errorCode = writeNode(bio, *node.getRightChild()))) goto read_saved_node_cleanup;
    }
    
read_saved_node_cleanup:
    
    if (errorCode != 0) {
      delete node.rightChild;
      delete node.leftChild;
      
      node.leftChild = NULL;
    }
    
    return errorCode;
  }
  
  int readTree(misc_binaryIO* bio, dbarts::Tree& tree, const dbarts::Data& data, const size_t* treeIndices)
  {
    tree.top.clear();
    return readNode(bio, tree.top, data, treeIndices);
  }
  
  int readTree(misc_binaryIO* bio, dbarts::SavedTree& tree)
  {
    tree.top.clear();
    return readNode(bio, tree.top);
  }
}
