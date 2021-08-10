#ifndef QuenchInfo_H
#define QuenchInfo_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {

  struct QuenchInfo {
    std::string quenchType;
    double perturbationStrength;
    InfoTargetState infoTargetState;

    std::string getInfoForFilename() const {
      std::stringstream info;
      info << "quenchType_" << quenchType;
      info << "_perturbationStrength_" << perturbationStrength << "_";
      info << infoTargetState.getInfoForFilename() << "/";
      return info.str();
    }

  };

}

#endif
