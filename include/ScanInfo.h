#ifndef ScanInfo_H
#define ScanInfo_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {

  struct ScanInfo {
    LiebLinState seedState;
    LiebLinState refStateAdditionalOutput;
    std::string selectionCriterion;
    std::string additionalOutput;
    double maxAbsWeight;
    std::pair<int, int> changeDoubledQuantumNumbs;
    double perturbationStrength;

    std::string getInfoForFilename() const {
      std::stringstream info;
      info << seedState.getInfoForFilename();
      info << "selectionCriterion_" << selectionCriterion;
      info << "_maxAbsWeight_" << maxAbsWeight;
      info << "_additionalOutput_" << additionalOutput;
      info << "_changeDoubledQuantumNumbs" << changeDoubledQuantumNumbs.first;
      info << "-" << changeDoubledQuantumNumbs.second;
      info << "/";
      mkdir(info.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
      return info.str();
    }
  };

}

#endif
