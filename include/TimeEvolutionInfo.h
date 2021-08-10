#ifndef TimeEvolutionInfo_H
#define TimeEvolutionInfo_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {

  class TimeEvolutionInfo {
    private:
      double totalTime;
      double timeStepSize;
      int timeSteps;
      std::string timeEvolvedOperator;

    public:
      int getTimeSteps() const {
        return timeSteps;
      };
      double getTimeStepSize() const {
        return timeStepSize;
      };
      double getTotalTime() const {
        return totalTime;
      };
      std::string getTimeEvolvedOperator() const {
        return timeEvolvedOperator;
      };
      std::string getInfoForFilename() const {
        std::stringstream info;
        info << "totalTime_" << totalTime;
        info << "_timeStepSize_" << timeStepSize << "_";
        return info.str();
      }

      TimeEvolutionInfo (const double& inputTotalTime, const double&
          inputTimeStepSize, const std::string& inputTimeEvolvedOperator)
        : totalTime(inputTotalTime),
          timeStepSize(inputTimeStepSize),
          timeEvolvedOperator(inputTimeEvolvedOperator) {
        timeSteps = std::floor(totalTime / timeStepSize);
      };

  };

}

#endif
