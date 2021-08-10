#include "ABACUS_VB.h"

namespace ABACUS_VB {

  NRG::NRG(const QuenchInfo& inputQuenchInfo, const ScanInfo& inputScanInfo,
      const NRGInfo& inputNRGInfo, const TimeEvolutionInfo&
      inputTimeEvolutionInfo)
    : quenchInfo(inputQuenchInfo),
      scanInfo(inputScanInfo),
      nrgInfo(inputNRGInfo),
      timeEvolutionInfo(inputTimeEvolutionInfo),
      currentNRGStep(0),
      currentBasisSize(nrgInfo.statesKept + nrgInfo.statesAdded) {
    (*this).setFilenameEnergyOutput();
    if (!(*this).hasFinishedBefore()) {
      (*this).loadComputationalBasis();
      (*this).initEnergiesComputationalBasis();
      (*this).initApproxEigenstateEnergies();
      (*this).initApproxEigenstateExpansions();
      (*this).initSelectedApproxEigenstates();
      (*this).initPerturbingOperator();
      (*this).initNRGHamiltonian();
      (*this).processNRGHamiltonian();
      for (int i = 0; i < (*this).getMaxNRGSteps(); i++) {
        currentNRGStep++;
        currentBasisSize += nrgInfo.statesAdded;
        (*this).updatePerturbingOperator();
        (*this).updateSelectedApproxEigenstates();
        (*this).updateNRGHamiltonian();
        (*this).processNRGHamiltonian();
      }
      (*this).setFilenameReturnAmplitudeAndFidelityOutput();
      (*this).setFilenameTimeEvolutionSelectedOperator();
      (*this).computeReturnAmplitudeAndFidelity();
      (*this).computeApproximateTimeEvolutionSelectedOperator();
    }
  }

}
