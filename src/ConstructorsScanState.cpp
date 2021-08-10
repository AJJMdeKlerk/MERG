#include "ABACUS_VB.h"

namespace ABACUS_VB {

  ScanState::ScanState(const ScanInfo& scanInfo)
    : doubledQuantumNumbsSeedState(scanInfo.seedState.getDoubledQuantumNumbs()),
      doubledQuantumNumbs(scanInfo.seedState.getDoubledQuantumNumbs()),
      foundMoverTypes(false),
      rightMoversPresent(false) {
    (*this).findMoverTypes();
    (*this).findIndexLeftmostRightmover();
    (*this).findIndexRightmostLeftmover();
  }

  ScanState::ScanState(const ScanInfo& scanInfo, const LiebLinState& state)
    : doubledQuantumNumbsSeedState(scanInfo.seedState.getDoubledQuantumNumbs()),
      doubledQuantumNumbs(state.getDoubledQuantumNumbs()),
      foundMoverTypes(false),
      rightMoversPresent(false) {
    (*this).findMoverTypes();
    (*this).findIndexLeftmostRightmover();
    (*this).findIndexRightmostLeftmover();
  }

  ScanState::ScanState(const ScanInfo& scanInfo, const std::string& label)
    : doubledQuantumNumbsSeedState(scanInfo.seedState.getDoubledQuantumNumbs()),
      doubledQuantumNumbs(convertAbsoluteLabelToDoubledQuantumNumbs(label, scanInfo)),
      foundMoverTypes(false),
      rightMoversPresent(false) {
    (*this).findMoverTypes();
    (*this).findIndexLeftmostRightmover();
    (*this).findIndexRightmostLeftmover();
  }

  ScanState::ScanState(const ScanInfo& scanInfo, const Eigen::VectorXi&
      inputDoubledQuantumNumbs)
    : doubledQuantumNumbsSeedState(scanInfo.seedState.getDoubledQuantumNumbs()),
      doubledQuantumNumbs(inputDoubledQuantumNumbs),
      foundMoverTypes(false),
      rightMoversPresent(false) {
      doubledQuantumNumbsSeedState =
        scanInfo.seedState.getDoubledQuantumNumbs();
      doubledQuantumNumbs = inputDoubledQuantumNumbs;
      (*this).findMoverTypes();
      (*this).findIndexLeftmostRightmover();
      (*this).findIndexRightmostLeftmover();
  }

  ScanState::ScanState() {}

}
