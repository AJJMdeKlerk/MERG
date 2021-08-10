#include "ABACUS_VB.h"

namespace ABACUS_VB {

  void ScanState::setFromScanInfo(const ScanInfo& scanInfo) {
    foundMoverTypes = false;
    rightMoversPresent = false;
    doubledQuantumNumbsSeedState = scanInfo.seedState.getDoubledQuantumNumbs();
    doubledQuantumNumbs = scanInfo.seedState.getDoubledQuantumNumbs();
    (*this).findMoverTypes();
    (*this).findIndexLeftmostRightmover();
    (*this).findIndexRightmostLeftmover();
  }

  std::string ScanState::getLabelAndWeight(const ScanInfo& scanInfo) const {
    std::stringstream labelAndWeight;
    std::cout.precision(6);
    labelAndWeight << (*this).getAbsoluteLabel() << ", " << (*this).getWeight(scanInfo); 
    return labelAndWeight.str();
  }

  double ScanState::getWeight(const ScanInfo& scanInfo) const {
    double weight;
    if (scanInfo.selectionCriterion == "Energy") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      weight = state.getEnergy();
    }
    else if (scanInfo.selectionCriterion == "EnergyDifference") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      weight = state.getEnergy() - scanInfo.seedState.getEnergy();
    }
    else if (scanInfo.selectionCriterion == "InteractionQuench") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      weight = std::abs(state.getEnergy() - scanInfo.seedState.getEnergy())
        / std::abs(scanInfo.perturbationStrength
          * G2(state, scanInfo.seedState));
    }
    else if (scanInfo.selectionCriterion == "InverseRho") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      weight = 1.0 / std::abs(rho(state, scanInfo.seedState));
    }
    else if (scanInfo.selectionCriterion == "QuantumNumber") {
      weight = std::max(doubledQuantumNumbs.maxCoeff(),
          std::abs(doubledQuantumNumbs.minCoeff()));
    }
    else {
      std::cerr << "Error, invalid selection criterion in getWeight. Exiting. ";
      std::cerr << std::endl;
      exit(-1);
    }
    return weight;
  }

  std::string ScanState::getAbsoluteLabel() const {
    std::stringstream binaryPartLabel;
    int particlesFound = 0;
    int position = doubledQuantumNumbs(0);
    while (particlesFound < doubledQuantumNumbs.size()) {
      if (position == doubledQuantumNumbs(particlesFound)) {
        binaryPartLabel << '1';
        particlesFound++;
      }
      else {
        binaryPartLabel << '0';
      }
      position += 2;
    }
    std::stringstream absoluteLabel;
    absoluteLabel << doubledQuantumNumbs(0) << ":";
    absoluteLabel << convertBinaryToHex(binaryPartLabel.str());
    return absoluteLabel.str();
  }

  std::string ScanState::getLabel() const {
    std::stringstream labelLeftmovers;
    std::stringstream labelRightmovers;
    for(int i = 0; i < doubledQuantumNumbs.size(); i++) {
      if(doubledQuantumNumbs(i) < doubledQuantumNumbsSeedState(i)) {
        if(!labelLeftmovers.str().empty())
          labelLeftmovers << ":";
        labelLeftmovers << doubledQuantumNumbsSeedState(i);
        labelLeftmovers << "@" << doubledQuantumNumbs(i);
      }
      else if(doubledQuantumNumbs(i) > doubledQuantumNumbsSeedState(i)) {
        if(!labelRightmovers.str().empty())
          labelRightmovers << ":";
        labelRightmovers << doubledQuantumNumbsSeedState(i);
        labelRightmovers << "@" << doubledQuantumNumbs(i);
      }
    }
    return labelLeftmovers.str() + "|" + labelRightmovers.str();
  }

  std::string ScanState::getAdditionalOutput(const ScanInfo& scanInfo) const {
    std::stringstream additionalOutput;
    additionalOutput << std::setprecision(16);
    if (scanInfo.additionalOutput == "Energy") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      additionalOutput << ", " << state.getEnergy();
    }
    else if (scanInfo.additionalOutput == "EnergyAndRho") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      additionalOutput << ", " << state.getEnergy();
      additionalOutput << ", " << rho(scanInfo.refStateAdditionalOutput, state);
    }
    else if (scanInfo.additionalOutput == "G2ME") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      additionalOutput << ", " << G2(scanInfo.refStateAdditionalOutput, state);
    }
    else if (scanInfo.additionalOutput == "Rho") {
      LiebLinState state(scanInfo.seedState.getsystemInfo(),
        doubledQuantumNumbs);
      additionalOutput << ", " << rho(scanInfo.refStateAdditionalOutput, state);
    }
    else if (scanInfo.additionalOutput == "SumIx2") {
      int sumDoubledQuantumNumbs = (*this).getSumDoubledQuantumNumbs();
      additionalOutput << ", " << sumDoubledQuantumNumbs;
    }
    else {
      std::cerr << "Error, invalid additional output criteria given. Exiting. ";
      exit(-1);
    }
    return additionalOutput.str();
  }

  int ScanState::getSumDoubledQuantumNumbs() const {
    return doubledQuantumNumbs.sum();
  }

  int ScanState::getNumberOfDoubledQuantumNumbs() const {
    return doubledQuantumNumbs.size();
  }

  Eigen::VectorXi ScanState::getDoubledQuantumNumbsSeedState() const {
    return doubledQuantumNumbsSeedState;
  }

  Eigen::VectorXi ScanState::getDoubledQuantumNumbs() const {
    return doubledQuantumNumbs;
  }

  Eigen::VectorXi ScanState::getMoverTypes() const {
    return moverTypes;
  }

}
