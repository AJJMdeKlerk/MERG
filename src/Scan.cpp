#include "ABACUS_VB.h"

namespace ABACUS_VB {

  std::string Scan::getFilenameApprovedOutput() const {
    return filenameApprovedOutput;
  }

  Scan::Scan(const ScanInfo& inputScanInfo)
    : scanInfo(inputScanInfo) {
    seedScanState.setFromScanInfo(scanInfo);
    (*this).setFilenameApprovedOutput();
    (*this).setFilenameDiscardedOutput();
    if (!(*this).isPerformedBefore()) {
      (*this).generateMomentumSectorSeedStates();

      std::cout << "The quantum numbers of the seed states for the different ";
      std::cout << "momentum sectors are: " << std::endl;
      for (auto momentumSectorSeedState : momentumSectorSeedStates) {
        std::cout << momentumSectorSeedState.getDoubledQuantumNumbs().transpose();
        std::cout << std::endl;
      }
      std::cout << std::endl;

      (*this).initApprovedOutput();
      (*this).generateOutput();
      (*this).writeOutputToDisc();
    }
  }

  void Scan::generateMomentumSectorSeedStates() {
    ScanState leftmovingSeedState(seedScanState);
    ScanState rightmovingSeedState(seedScanState);
    int numberOfSeedStatesGenerated = 0;
    if (seedScanState.isApprovedByWeightAndMomentum(scanInfo)) {
      momentumSectorSeedStates.push_back(seedScanState);
      numberOfSeedStatesGenerated++;
    }
    int indexLastParticle = seedScanState.getNumberOfDoubledQuantumNumbs() - 1;
    int newStatesGenerated = -1;
    int i = 0;
    while (numberOfSeedStatesGenerated < 1 || newStatesGenerated != 0) {
      newStatesGenerated = 0;
      leftmovingSeedState.leftmoveDoubledQuantumNumber(i);
      leftmovingSeedState.resetMoverInfo();
      rightmovingSeedState.rightmoveDoubledQuantumNumber(indexLastParticle - i);
      rightmovingSeedState.resetMoverInfo();
      if (leftmovingSeedState.isApprovedByWeightAndMomentum(scanInfo)) {
        momentumSectorSeedStates.push_back(leftmovingSeedState);
        numberOfSeedStatesGenerated++;
        newStatesGenerated++;
      }
      if (rightmovingSeedState.isApprovedByWeightAndMomentum(scanInfo)) {
        momentumSectorSeedStates.push_back(rightmovingSeedState);
        numberOfSeedStatesGenerated++;
        newStatesGenerated++;
      }
      i = (i + 1) % rightmovingSeedState.getNumberOfDoubledQuantumNumbs();
    }
  }

  void Scan::generateOutput() {
    (*this).generateBetheTreeWithCutoff(momentumSectorSeedStates);
    (*this).startRoundsOfForcedDescents();
  }

  void Scan::generateBetheTreeWithCutoff(std::vector<ScanState> seeds) {
    tbb::parallel_for_each(seeds,
      [&](std::vector<ScanState>::reference state,
        tbb::feeder<ScanState>& feeder) {
        std::queue<ScanState> descendents =
          state.generateEqualMomentumDescendents();
        while(!descendents.empty()) {
          std::string label = descendents.front().getAbsoluteLabel();
          double weight = descendents.front().getWeight(scanInfo);
          if (isApprovedByWeight(scanInfo, weight)) {
            approvedOutput.push(label + ", " + std::to_string(weight) +
              descendents.front().getAdditionalOutput(scanInfo));
            feeder.add(descendents.front());
          }
          else
            pausedBranches.push(label);
          descendents.pop();
        }
      }
    );
  }

  void Scan::startRoundsOfForcedDescents() {
    bool keepDescending = true;
    int numberOfForcedDescents = 1;
    while (keepDescending) {
      // unsafe_size is actually safe in this scenario since we are done
      // accessing the concurrent container it acts on when we invoke it
      auto approvedStatesBeforeDescent = approvedOutput.unsafe_size() - 1;
      (*this).forceDescent();
      auto approvedStatesAfterDescent = approvedOutput.unsafe_size() - 1;
      std::cout << "Round " << numberOfForcedDescents;
      std::cout << " of forced descents generated ";
      std::cout << approvedStatesAfterDescent - approvedStatesBeforeDescent;
      std::cout << " new approved states. " << std::endl;
      std::cout << "This makes a total of " << approvedStatesAfterDescent;
      std::cout << " approved states thus far. " << std::endl;
      std::cout << "The next round, if any, will start with ";
      std::cout << pausedBranches.unsafe_size() << " seed states.";
      std::cout << std::endl << std::endl;
      if (approvedStatesBeforeDescent == approvedStatesAfterDescent)
        keepDescending = false;
      else
        numberOfForcedDescents++;
    }
    std::cout << "Round " << numberOfForcedDescents;
    std::cout << " of forced descents generated no new approved states. ";
    std::cout << "This ends this part of the procedure. ";
    std::cout << std::endl << std::endl;
  }

  void Scan::forceDescent(){
    auto labelsSeedsForcedDescent = (*this).getLabelsSeedsForcedDescent();
    tbb::parallel_for_each(labelsSeedsForcedDescent, 
      [&](std::vector<std::string>::reference labelParent,
        tbb::feeder<std::string>& feeder) {
        ScanState parent(scanInfo, labelParent);
        double weightParent = parent.getWeight(scanInfo);
        auto descendents = parent.generateEqualMomentumDescendents();
        while(!descendents.empty()) {
          auto labelDescendent = descendents.front().getAbsoluteLabel();
          auto weightDescendent = descendents.front().getWeight(scanInfo);
          if (isApprovedByWeight(scanInfo, weightDescendent)) {
            approvedOutput.push(labelDescendent + ", " +
                std::to_string(weightDescendent) +
                descendents.front().getAdditionalOutput(scanInfo));
            feeder.add(labelDescendent);
          }
          else if (std::abs(weightDescendent) < std::abs(weightParent)) {
            feeder.add(labelDescendent);
          }
          else
            pausedBranches.push(labelDescendent);
          descendents.pop();
        }
      }
    );
  }

  std::vector<std::string> Scan::getLabelsSeedsForcedDescent() {
    std::vector<std::string> labelsSeedsForcedDescent;
    while (!pausedBranches.empty()) {
      std::string labelPausedState;
      if (pausedBranches.try_pop(labelPausedState))
        labelsSeedsForcedDescent.push_back(labelPausedState);
    }
    return labelsSeedsForcedDescent;
  }

  bool Scan::isPerformedBefore() {
    std::ifstream outputFile;
    outputFile.open(filenameApprovedOutput);
    if (outputFile.peek() != std::ifstream::traits_type::eof()) {
      std::cout << "The file named " << filenameApprovedOutput << " already ";
      std::cout << "existed, so this scanning routine has already been ";
      std::cout << "performed and we will not repeat it. " << std::endl;
      std::cout << std::endl;
      return true;
    }
    else
      return false;
  }

  void Scan::writeOutputToDisc() {
    writeStringsToDisc(approvedOutput, filenameApprovedOutput);
  }

  void Scan::initApprovedOutput() {
    std::string columns = (*this).getColumnsApprovedOutput();
    approvedOutput.push(columns);
    for(unsigned int i = 0; i < momentumSectorSeedStates.size(); i++) {
      approvedOutput.push(momentumSectorSeedStates[i].getLabelAndWeight(scanInfo)
        + momentumSectorSeedStates[i].getAdditionalOutput(scanInfo));
    }
  }

  std::string Scan::getColumnsApprovedOutput() {
    std::string columns = (*this).getColumnsDiscardedOutput();
    if (scanInfo.additionalOutput == "")
      return columns;
    else if (scanInfo.additionalOutput == "Energy")
      return columns + ", Energy";
    else if (scanInfo.additionalOutput == "EnergyAndRho")
      return columns + ", Energy, Rho";
    else if (scanInfo.additionalOutput == "G2ME")
      return columns + ", G2ME";
    else if (scanInfo.additionalOutput == "Rho")
      return columns + ", Rho";
    else if (scanInfo.additionalOutput == "SumIx2")
      return columns + ", SumIx2";
    else {
      std::cerr << "Error, invalid additional output criteria given. Exiting. ";
      std::cerr << std::endl;
      exit(-1);
    }
  }

  std::string Scan::getColumnsDiscardedOutput() {
    return "Label, " + scanInfo.selectionCriterion;
  }

  void Scan::setFilenameApprovedOutput() {
    std::stringstream name;
    mkdir(scanInfo.getInfoForFilename().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << scanInfo.getInfoForFilename();
    name << "ApprovedStates.csv";
    filenameApprovedOutput = name.str();
  }

  void Scan::setFilenameDiscardedOutput() {
    std::stringstream name;
    mkdir(scanInfo.getInfoForFilename().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << scanInfo.getInfoForFilename();
    name << "DiscardedStates.csv";
    filenameDiscardedOutput = name.str();
  }

  void writeStringsToDisc(tbb::concurrent_queue<std::string> strings,
      const std::string& filename) {
    std::string singleString;
    std::ofstream outputFile(filename);
    while (!strings.empty()) {
      strings.try_pop(singleString);
      outputFile << singleString << '\n';
    }
    outputFile.close();
  }

  std::vector<std::pair<std::string, double> > readWeightedLabelsFromDisc
    (const std::string& filename) {
    std::ifstream input(filename);
    std::vector<std::pair<std::string, double> > weightedLabels;
    std::string columns;
    std::string label;
    std::string remainingOutput;
    double weight;
    getline(input, columns);
    while (getline(input, label, ',')) {
      getline(input, remainingOutput);
      std::stringstream splittableRemainingOutput(remainingOutput);
      splittableRemainingOutput >> weight;
      weightedLabels.push_back(std::make_pair(label, weight));
    }
    input.close();
    return weightedLabels;
  }

  std::vector<LiebLinState> setStatesFromWeightedLabels
    (const std::vector<std::pair<std::string, double> >& weightedLabels,
     const SystemInfo& systemInfo, const ScanInfo& scanInfo){
    int numberOfStates = weightedLabels.size();
    std::vector<LiebLinState> setStates;
    setStates.resize(numberOfStates);
    tbb::parallel_for (0, numberOfStates,
      [&](int i) {
        LiebLinState liebLinState(systemInfo, scanInfo,
          weightedLabels[i].first);
        setStates[i] = liebLinState;
      }
    );
    return setStates;
  }

}
