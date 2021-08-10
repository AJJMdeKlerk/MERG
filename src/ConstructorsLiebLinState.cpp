#include "ABACUS_VB.h"

namespace ABACUS_VB {

  LiebLinState::LiebLinState() {}

  LiebLinState::LiebLinState(const SystemInfo& inputScanInfo, const
      ScanInfo& ScanInfo, const std::string& label)
  : systemInfo(inputScanInfo),
    allComputed(false) {
    doubledQuantumNumbs = convertAbsoluteLabelToDoubledQuantumNumbs(label, ScanInfo);
    (*this).defaultInitAccuracy();
    (*this).checkPositiveInteractionStrength();
    (*this).defaultInitRapidities();
    (*this).computeAll();
  }

  LiebLinState::LiebLinState (const SystemInfo& inputsystemInfo)
  : systemInfo(inputsystemInfo),
    allComputed(false) {
    (*this).defaultInitAccuracy();
    (*this).checkPositiveInteractionStrength();
    (*this).defaultInitDoubledQuantumNumbs();;
    (*this).defaultInitRapidities();
    (*this).computeAll();
  }

  LiebLinState::LiebLinState (const SystemInfo& inputsystemInfo,
    const Eigen::VectorXi& inputDoubledQuantumNumbers)
  : systemInfo(inputsystemInfo),
    allComputed(false) {
    doubledQuantumNumbs = inputDoubledQuantumNumbers;
    (*this).defaultInitAccuracy();
    (*this).checkPositiveInteractionStrength();
    (*this).defaultInitRapidities();
    (*this).computeAll();
  }

  void LiebLinState::defaultInitDoubledQuantumNumbs() {
    // By default we initialise in the ground state
  	doubledQuantumNumbs.resize((*this).getParticleNumber());
    for (int i = 0; i < (*this).getParticleNumber(); i++) {
      doubledQuantumNumbs(i) = - ((*this).getParticleNumber() + 1);
      doubledQuantumNumbs(i) += 2 * (i + 1);
    }
  }

  void LiebLinState::defaultInitRapidities() {
  	rapidities.resize((*this).getParticleNumber());
  	for (int i = 0; i < (*this).getParticleNumber(); i++) {
  		rapidities(i) = M_PI * doubledQuantumNumbs(i)
        / (*this).getSystemLength();
  	}
  }

  void LiebLinState::defaultInitAccuracy() {
    maxErrorRapidities = (*this).getParticleNumber() * 1.0e-10;
    maxNewtonIterations = 100;
  }

}
