#ifndef Scan_H
#define Scan_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {

  class Scan {
    private:
      ScanInfo scanInfo;
      ScanState seedScanState;

      std::vector<ScanState> momentumSectorSeedStates;
      tbb::concurrent_queue<std::string> approvedOutput;
      tbb::concurrent_queue<std::string> pausedBranches;
      std::string filenameApprovedOutput;
      std::string filenameDiscardedOutput;
      void setFilenameApprovedOutput();
      void setFilenameDiscardedOutput(); // obsolete
      std::string getColumnsApprovedOutput();
      std::string getColumnsDiscardedOutput(); // obsolete

      bool isPerformedBefore();
      void initApprovedOutput();
      void generateOutput(); 
      void generateMomentumSectorSeedStates(); 
      void generateBetheTreeWithCutoff(std::vector<ScanState> seeds);
      void startRoundsOfForcedDescents();
      void forceDescent();
      void writeOutputToDisc();
      std::vector<std::string> getLabelsSeedsForcedDescent();
    public:
      Scan(const ScanInfo& inputParameters);
      std::string getFilenameApprovedOutput() const;
  };

  void writeStringsToDisc(tbb::concurrent_queue<std::string> strings, 
      const std::string& filename);
  std::vector<std::pair<std::string, double> > readWeightedLabelsFromDisc
    (const std::string& filename);
  std::vector<LiebLinState> setStatesFromWeightedLabels 
    (const std::vector<std::pair<std::string, double> >& weightedLabels, 
     const SystemInfo& systemInfo, const ScanInfo& scanInfo);

}

#endif
