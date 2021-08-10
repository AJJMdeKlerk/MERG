#include "ABACUS_VB.h"

namespace ABACUS_VB {

  std::string NRG::getFilenameEnergyOutput() const {
    return filenameEnergyOutput;
  }

  std::string NRG::getFilenameReturnAmplitudeAndFidelityOutput() const {
    return filenameReturnAmplitudeAndFidelityOutput;
  }

  std::string NRG::getFilenameTimeEvolutionSelectedOperator() const {
    return filenameTimeEvolutionSelectedOperator;
  }

  int NRG::getMaxNRGSteps() const {
    return nrgInfo.maxNRGSteps;
  }

  int NRG::getCurrentNRGStep() const {
    return currentNRGStep;
  }

  int NRG::getCurrentBasisSize() const {
    return nrgInfo.statesKept + nrgInfo.statesAdded 
      + (*this).getCurrentNRGStep() * nrgInfo.statesAdded;
  }

  int NRG::getMaxBasisSize() const {
    return nrgInfo.statesKept + nrgInfo.statesAdded 
      + (*this).getMaxNRGSteps() * nrgInfo.statesAdded;
  }

  double NRG::getTargetEnergy() const {
    return quenchInfo.infoTargetState.energy;
  }

}
