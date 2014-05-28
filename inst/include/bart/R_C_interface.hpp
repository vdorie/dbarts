#ifndef BART_R_C_INTERFACE_H
#define BART_R_C_INTERFACE_H

namespace bart {
  struct Control;
  struct Model;
  struct Data;
  
  struct CGMPrior;
  struct NormalPrior;
  struct ChiSquaredPrior;
  
  struct BARTFit;
  struct Results;
}

extern "C" {
  // pair calls of create<->destroy, initialize<->invalidate
  bart::BARTFit* bart_createFit(bart::Control* control, bart::Model* model, bart::Data* data);
  void bart_initializeFit(bart::BARTFit* fit, bart::Control* control, bart::Model* model, bart::Data* data);
  void bart_destroyFit(bart::BARTFit* fit);
  void bart_invalidateFit(bart::BARTFit* fit);
  
  bart::Results* bart_runSampler(bart::BARTFit* fit);
  void bart_setResponse(bart::BARTFit* fit, const double* newResponse);
  
  
  bart::CGMPrior* bart_createCGMPrior();
  bart::CGMPrior* bart_createCGMPriorFromControl(const bart::Control* control);
  void bart_destroyCGMPrior(bart::CGMPrior* prior);
  void bart_initializeCGMPriorFromControl(bart::CGMPrior* prior, const bart::Control* control);
  void bart_invalidateCGMPrior(bart::CGMPrior* prior);
  
  bart::NormalPrior* bart_createNormalPrior();
  bart::NormalPrior* bart_createNormalPriorFromControl(const bart::Control* control);
  void bart_destroyNormalPrior(bart::NormalPrior* prior);
  void bart_initializeNormalPriorFromControl(bart::NormalPrior* prior, const bart::Control* control);
  void bart_invalidateNormalPrior(bart::NormalPrior* prior);
  
  bart::ChiSquaredPrior* bart_createChiSquaredPrior();
  bart::ChiSquaredPrior* bart_createChiSquaredPriorFromControl(const bart::Control* control);
  void bart_destroyChiSquaredPrior(bart::ChiSquaredPrior* prior);
  void bart_initializeChiSquaredPriorFromControl(bart::ChiSquaredPrior* prior, const bart::Control* control);
  void bart_invalidateChiSquaredPrior(bart::ChiSquaredPrior* prior);
}
#endif // BART_R_C_INTERFACE_H
