#ifndef NRGInfo_H
#define NRGInfo_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {
  
  struct NRGInfo {
    int statesKept;
    int statesAdded;
    int maxNRGSteps;
    std::string method;

    int getTotalStates() const {
      return statesKept + (maxNRGSteps + 1) * statesAdded;
    }

    std::string getInfoForFilename() const {
      std::stringstream info;
      info << "statesKept" << statesKept;
      info << "_statesAdded_" << statesAdded;
      info << "_maxNRGSteps_" << maxNRGSteps;
      info << "_method" << method << "/";
      return info.str();
    }

  };

}

#endif
