#ifndef InfoTargetState_H
#define InfoTargetState_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {
  
  struct InfoTargetState {
    double energy;

    std::string getInfoForFilename() const {
      std::stringstream info;
      info << "energyTargetState_" << energy;
      return info.str();
    }

  };

}

#endif
