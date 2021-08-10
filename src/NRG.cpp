#include "ABACUS_VB.h"

namespace ABACUS_VB {

  void NRG::processNRGHamiltonian() {
    (*this).diagonaliseHamiltonianAndUpdateApproxEigenvectors();
    (*this).writeResultsDiagonalisationToTerminal();
    (*this).writeResultsDiagonalisationToFile();
  }
 
  void NRG::diagonaliseHamiltonianAndUpdateApproxEigenvectors() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ces(nrgHamiltonian);
    selectedApproxEigenstateExpansions *= ces.eigenvectors().real();
    std::vector<std::pair<double, int> > indexedInnerProducts(nrgInfo.statesKept
        + nrgInfo.statesAdded); 
    tbb::parallel_for (0, nrgInfo.statesKept + nrgInfo.statesAdded,
      [&](int i) {
        double innerProduct;
        innerProduct = selectedApproxEigenstateExpansions.col(i).transpose() *
          approxEigenstateExpansions.col(0).head(currentBasisSize);
        indexedInnerProducts[i].first = 1.0 / std::abs(innerProduct);
        indexedInnerProducts[i].second = i;
      }
    );
    std::sort(indexedInnerProducts.begin(), indexedInnerProducts.end());
    if (1.0 / indexedInnerProducts[0].first < 0.75) {
      std::cerr << "Largest overlap between NRG-steps smaller than 0.75. ";
      std::cerr << "Exiting to prevent possibly following the wrong ";
      std::cerr << "approximate eigenstate. Change to smaller NRG-steps.";
      exit(-1);
    }
    tbb::parallel_for(0, nrgInfo.statesKept + nrgInfo.statesAdded,
      [&](int i) {
        approxEigenstateExpansions.col(selectedEigenstates(i)).head(currentBasisSize)
          = selectedApproxEigenstateExpansions.col(indexedInnerProducts[i].second);
        approxEigenstateEnergies(selectedEigenstates(i)) 
          = ces.eigenvalues()(indexedInnerProducts[i].second);
      }
    );
  }

  void NRG::writeResultsDiagonalisationToTerminal() {
    std::cout << "The targeted eigenvalue of H_NRG after NRGstep ";
    std::cout << std::setprecision(16) << currentNRGStep;
    std::cout << ", which corresponds to having taken into account ";
    std::cout << currentBasisSize;
    std::cout << " states is: " << approxEigenstateEnergies(0) << std::endl;
  }

  void NRG::writeResultsDiagonalisationToFile() {
    std::ofstream statesVsEnergy(filenameEnergyOutput, std::fstream::app);
    statesVsEnergy << std::setprecision(16);
    statesVsEnergy << currentBasisSize << ", ";
    statesVsEnergy << approxEigenstateEnergies(0) / (*this).getTargetEnergy() - 1;
    statesVsEnergy << ", " << 1.0 / weightedLabels[currentBasisSize].second;
    statesVsEnergy << std::endl;
    statesVsEnergy.close();
  }

  void NRG::updateNRGHamiltonian() {
    nrgHamiltonian.noalias() = selectedApproxEigenstateExpansions.transpose() *
      perturbingOperator.topLeftCorner(currentBasisSize, currentBasisSize) *
      selectedApproxEigenstateExpansions;
    nrgHamiltonian.noalias() += selectedApproxEigenstateExpansions.transpose() *
      energiesComputationalBasis.head(currentBasisSize).asDiagonal()
      * selectedApproxEigenstateExpansions;
  }

  void NRG::updateSelectedApproxEigenstates() {
    (*this).computeWeightsAndSelectApproxEigenstates();
    (*this).addNewStatesToSelectedApproxEigenstates();
    (*this).computeSelectedApproxEigenstateExpansions();
  }

  void NRG::computeWeightsAndSelectApproxEigenstates() {
    std::vector<std::pair<double, int> > weightedApproxEigenstates;
    if (nrgInfo.method == "TSA" || nrgInfo.method == "HOSTS") {
      weightedApproxEigenstates.resize(nrgInfo.statesKept + nrgInfo.statesAdded);
      for (int i = 1; i < nrgInfo.statesKept + nrgInfo.statesAdded; i++) { 
        weightedApproxEigenstates[i].first = 
          approxEigenstateEnergies(selectedEigenstates(i));
        weightedApproxEigenstates[i].second = selectedEigenstates(i);
      }
    }
    else if (nrgInfo.method == "MERG") {
      Eigen::MatrixXd
        MEsPerturbingOperatorNewAndApproxEigenstates(nrgInfo.statesAdded,
            currentBasisSize - nrgInfo.statesAdded);
      MEsPerturbingOperatorNewAndApproxEigenstates.noalias() =
        perturbingOperator.block(currentBasisSize - nrgInfo.statesAdded,  0,
            nrgInfo.statesAdded, currentBasisSize - nrgInfo.statesAdded)
        * approxEigenstateExpansions.topLeftCorner(currentBasisSize -
            nrgInfo.statesAdded, currentBasisSize - nrgInfo.statesAdded);
      weightedApproxEigenstates.resize(currentBasisSize - nrgInfo.statesAdded);
      tbb::parallel_for(1, currentBasisSize - nrgInfo.statesAdded,
          [&](int i) {
          double weight = 0.0;
          for (int j = 0; j < nrgInfo.statesAdded; j++) {
            weight += (approxEigenstateEnergies[0] - approxEigenstateEnergies[i])
              * (approxEigenstateEnergies[0] - approxEigenstateEnergies[j])
              / (MEsPerturbingOperatorNewAndApproxEigenstates(j, 0) 
                  * MEsPerturbingOperatorNewAndApproxEigenstates(j, i));
          }
          weightedApproxEigenstates[i].first = std::abs(weight);
          weightedApproxEigenstates[i].second = i;
        }
      );
    }
    else {
      std::cerr << "The NRG-method you requested has currently not been ";
      std::cerr << "implemented yet. Exiting";
      std::cerr << std::endl;
      exit(-1);
    }
    weightedApproxEigenstates[0].first = - 1.0;
    weightedApproxEigenstates[0].second = 0;
    std::sort(weightedApproxEigenstates.begin(),
        weightedApproxEigenstates.end());
    for(int i = 0; i < nrgInfo.statesKept; i++)
      selectedEigenstates(i) = weightedApproxEigenstates[i].second; 
  }

  void NRG::addNewStatesToSelectedApproxEigenstates() {
    for (int i = nrgInfo.statesKept; i < nrgInfo.statesKept +
        nrgInfo.statesAdded; i++) { 
      selectedEigenstates(i) = currentBasisSize -
        nrgInfo.statesKept - nrgInfo.statesAdded + i; 
    }
  }
  
  void NRG::computeSelectedApproxEigenstateExpansions() {
    selectedApproxEigenstateExpansions.resize(currentBasisSize,
        nrgInfo.statesKept + nrgInfo.statesAdded);
    for (int i = 0; i < nrgInfo.statesKept + nrgInfo.statesAdded; i++) {
      selectedApproxEigenstateExpansions.col(i) =
        approxEigenstateExpansions.col(selectedEigenstates(i)).head(currentBasisSize);
    }
  }

  void NRG::updatePerturbingOperator() {
    tbb::parallel_for(currentBasisSize - nrgInfo.statesAdded, currentBasisSize,
      [&](int j) {
        for (int i = 0; i <= j; i++) {
          double matrixElement = (*this).getMEPerturbingOperator(i, j);
          perturbingOperator(i, j) = matrixElement;
          perturbingOperator(j, i) = matrixElement;
        }
      }
    );
  }

  void NRG::initApproxEigenstateExpansions() {
    approxEigenstateExpansions =
      Eigen::MatrixXd::Identity((*this).getMaxBasisSize(),
          (*this).getMaxBasisSize());
  }

  void NRG::initSelectedApproxEigenstates() {
    selectedEigenstates.resize(currentBasisSize);
    for (int i = 0; i < currentBasisSize; i++) {
      selectedEigenstates(i) = i;
    }
    (*this).computeSelectedApproxEigenstateExpansions();
  }

  void NRG::initNRGHamiltonian() {
    nrgHamiltonian = perturbingOperator.topLeftCorner(currentBasisSize,
        currentBasisSize); 
    nrgHamiltonian +=
      energiesComputationalBasis.head(currentBasisSize).asDiagonal();
  }

  void NRG::initPerturbingOperator() {
    perturbingOperator.resize((*this).getMaxBasisSize(),
        (*this).getMaxBasisSize());
    tbb::parallel_for(0, currentBasisSize,
      [&](int j) {
        for (int i = 0; i <= j; i++) {
          double matrixElement = (*this).getMEPerturbingOperator(i, j);
          perturbingOperator(i, j) = matrixElement;
          perturbingOperator(j, i) = matrixElement;
        }
      }
    );
  }

  void NRG::setFilenameEnergyOutput() {
    std::stringstream name;
    name << scanInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << quenchInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << nrgInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << "Energies.csv";
    filenameEnergyOutput = name.str();
  }

  void NRG::setFilenameReturnAmplitudeAndFidelityOutput() {
    std::stringstream name;
    name << scanInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << quenchInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << nrgInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << "ReturnAmplitudeAndFidelity_";
    name << "currentBasisSize_";
    name << currentBasisSize;
    name << ".csv";
    filenameReturnAmplitudeAndFidelityOutput = name.str();
  }

  void NRG::setFilenameTimeEvolutionSelectedOperator() {
    std::stringstream name;
    name << scanInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << quenchInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << nrgInfo.getInfoForFilename();
    mkdir(name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    name << "TimeEvolution_";
    name << "currentBasisSize_";
    name << currentBasisSize;
    name << ".csv";
    filenameTimeEvolutionSelectedOperator = name.str();
  }

  Eigen::MatrixXd NRG::getWeightedMatrixElementsTimeEvolution() const {
    Eigen::MatrixXd weightedMatrixElementsTimeEvolution;
    if (quenchInfo.quenchType == "InteractionStrength" 
        && timeEvolutionInfo.getTimeEvolvedOperator() == "G2") {
      weightedMatrixElementsTimeEvolution = (perturbingOperator /
          quenchInfo.perturbationStrength).cwiseProduct
        (approxEigenstateExpansions.col(0) *
         approxEigenstateExpansions.col(0).transpose());
    }
    else {
      std::cerr << "The chosen operator to time evolve is not implemented ";
      std::cerr << "for the quench at hand. " << std::endl;
      exit(-1);
    }
    return weightedMatrixElementsTimeEvolution;
  }

  std::vector<std::vector<int> > NRG::getIndicesSelectedSummandsTimeEvolution
      (const Eigen::MatrixXd& weightedMatrixElementsTimeEvolution) const {
    std::vector<std::vector<int> > result(currentBasisSize);
    for (int i = 0; i < currentBasisSize; i++) {
      double sumAbsValuesCol =
        weightedMatrixElementsTimeEvolution.col(i).cwiseAbs().sum();
      std::vector<std::pair<double, int> > 
        weightedIndicesColumn(currentBasisSize);
      for (int j = 0; j < currentBasisSize; j++) {
        weightedIndicesColumn[j] = std::make_pair(1.0 /
            std::abs(weightedMatrixElementsTimeEvolution(j, i)), j);
      }
      std::sort(weightedIndicesColumn.begin(), weightedIndicesColumn.end());
      double intermediateSumAbsValuesCol = 0.0;
      int k = 0;
      while (intermediateSumAbsValuesCol < 0.99 * sumAbsValuesCol) {
        intermediateSumAbsValuesCol += 1.0 / weightedIndicesColumn[k].first;
        k++;
      }
      std::vector<int> selectedIndicesColumn(k);
      for (int l = 0; l < k; l++) {
        selectedIndicesColumn[l] = weightedIndicesColumn[l].second;
      }
      result[i] = selectedIndicesColumn;
    }
    return result;
  }

  void NRG::computeApproximateTimeEvolutionSelectedOperator() {
    auto weightedMatrixElements =
      (*this).getWeightedMatrixElementsTimeEvolution();
    auto selectedSummands =
      (*this).getIndicesSelectedSummandsTimeEvolution(weightedMatrixElements);

    // We compute the sum of the amplitudes of the different states comprising 
    // the superposition of the state we constructed weighted by the relevant
    // phases for each step of the time evolution
    std::vector<std::complex<double> >
      timeEvolution(timeEvolutionInfo.getTimeSteps());
    tbb::parallel_for(0, timeEvolutionInfo.getTimeSteps(),
      [&](int i) {
        std::complex<double> sum = 0.0;
        for (int j = 0; j < currentBasisSize; j++){
          for (int k = 0; k < int(selectedSummands[j].size()); k++){
            sum += std::exp(- II * (energiesComputationalBasis[j] -
                  energiesComputationalBasis[selectedSummands[j][k]]) * double(i) *
                timeEvolutionInfo.getTimeStepSize()) 
              * weightedMatrixElements(selectedSummands[j][k], j);
          }
        }
        timeEvolution[i] = sum;
      }
    );
    
    std::ofstream timeEvolutionData(filenameTimeEvolutionSelectedOperator,
        std::fstream::app);
    for (int i = 0; i < timeEvolutionInfo.getTimeSteps(); i++){
      timeEvolutionData << i * timeEvolutionInfo.getTimeStepSize() << ", ";  
      timeEvolutionData << std::real(timeEvolution[i]) << std::endl;
    }
    timeEvolutionData.close();

  }

  void NRG::computeReturnAmplitudeAndFidelity() {
    // We compute the sum of the amplitudes of the different states comprising 
    // of the superposition of the state we constructed weighted by the relevant
    // phases for each step of the time evolution
    std::vector<std::complex<double> >
      sumOverlapsWithPhases(timeEvolutionInfo.getTimeSteps());
    tbb::parallel_for(0, timeEvolutionInfo.getTimeSteps(),
      [&](int i) {
        std::complex<double> sum = 0.0;
        for (int j = 0; j < currentBasisSize; j++){
          sum += std::pow(std::abs(approxEigenstateExpansions(j, 0)), 2.0) 
            * std::exp(- II * energiesComputationalBasis[j] * double(i)
              * timeEvolutionInfo.getTimeStepSize());
        }
        sumOverlapsWithPhases[i] = sum;
      }
    );
    std::ofstream returnAmplitudeAndFidelityOutput
      (filenameReturnAmplitudeAndFidelityOutput, std::fstream::app); 
    returnAmplitudeAndFidelityOutput << std::setprecision(16);
    for (int i = 0; i < timeEvolutionInfo.getTimeSteps(); i++){
      returnAmplitudeAndFidelityOutput << i * timeEvolutionInfo.getTimeStepSize();
      returnAmplitudeAndFidelityOutput << ", ";
      returnAmplitudeAndFidelityOutput << std::real(sumOverlapsWithPhases[i]);
      returnAmplitudeAndFidelityOutput << ", ";
      returnAmplitudeAndFidelityOutput << std::imag(sumOverlapsWithPhases[i]);
      returnAmplitudeAndFidelityOutput << ", ";
      returnAmplitudeAndFidelityOutput <<
        std::pow(std::abs(sumOverlapsWithPhases[i]), 2.0) << std::endl;
    }
    returnAmplitudeAndFidelityOutput.close();
  }

  bool NRG::hasFinishedBefore() {
    bool donePreviously = false;
    std::ifstream fileEnergyOutput((*this).getFilenameEnergyOutput());
    std::string outputLine;
    int numberOfLines = 0;
    while (getline(fileEnergyOutput, outputLine))
      numberOfLines++;
    fileEnergyOutput.close();
    if (numberOfLines == (*this).getMaxNRGSteps() + 2) {
      std::cout << "This run has already been completed and stored in: ";
      std::cout << (*this).getFilenameEnergyOutput() << std::endl << std::endl;
      donePreviously = true;
    }
    else {
      std::ofstream fileEnergyOutput(filenameEnergyOutput);
      fileEnergyOutput << "N, PercentualEnergyError, LowestWeightIncluded";
      fileEnergyOutput << std::endl;
      fileEnergyOutput.close();
      std::ofstream returnAmplitudeAndFidelityOutput
        (filenameReturnAmplitudeAndFidelityOutput);
      returnAmplitudeAndFidelityOutput << "Time, RealPartReturnAmplitude, ";
      returnAmplitudeAndFidelityOutput << "ImaginaryPartReturnAmplitude, ";
      returnAmplitudeAndFidelityOutput << "Fidelity";
      returnAmplitudeAndFidelityOutput << std::endl;
      returnAmplitudeAndFidelityOutput.close();
      std::ofstream timeEvolutionSelectedOperator
        (filenameTimeEvolutionSelectedOperator);
      timeEvolutionSelectedOperator << "Time, "; 
      timeEvolutionSelectedOperator << "ExpectationValueOperator";
      timeEvolutionSelectedOperator << std::endl;
      timeEvolutionSelectedOperator.close();
    }
    return donePreviously;
  }

  bool sortByWeight(const std::pair<std::string, double> &a, 
      const std::pair<std::string, double> &b) { 
    return (a.second < b.second); 
  } 

  void NRG::initApproxEigenstateEnergies() {
    approxEigenstateEnergies.resize(nrgInfo.getTotalStates());
    for (int i = 0; i < nrgInfo.getTotalStates(); i++)
      approxEigenstateEnergies(i) = computationalBasis[i].getEnergy();
  }

  void NRG::initEnergiesComputationalBasis() {
    energiesComputationalBasis.resize((*this).getMaxBasisSize());
    for (int i = 0; i < nrgInfo.getTotalStates(); i++)
      energiesComputationalBasis(i) = computationalBasis[i].getEnergy();
  }

  void NRG::loadComputationalBasis() {
    Scan scan(scanInfo);
    weightedLabels = readWeightedLabelsFromDisc(scan.getFilenameApprovedOutput());
    std::sort(weightedLabels.begin(), weightedLabels.end(), sortByWeight);
    (*this).checkSizeWeightedLabels();
    (*this).computeSizesSymmetryBlocks();
    (*this).reorderSymmetryBlocks();
    computationalBasis = setStatesFromWeightedLabels(weightedLabels,
        scanInfo.seedState.getSystemInfo(), scanInfo);
  }

  void NRG::checkSizeWeightedLabels() {
    if (int(weightedLabels.size()) < nrgInfo.getTotalStates()) {
      std::cerr << "Error, the loaded computational basis is too small to ";
      std::cerr << "complete the NRG-procedure. Please change either the scan ";
      std::cerr << "or the number of NRG-steps." << std::endl;
      exit(-1);
    }
  }

  double NRG::getMEPerturbingOperator(const int& i, const int& j) const {
    if (quenchInfo.quenchType == "InteractionStrength") {
      return quenchInfo.perturbationStrength * 
        G2(computationalBasis[i], computationalBasis[j]);
    }
    else {
      std::cerr << "The quench you requested has currently not been ";
      std::cerr << "implemented yet. Your only option is: InteractionStrength";
      std::cerr << std::endl;
      exit(-1);
    }
  }

  void NRG::computeSizesSymmetryBlocks() {
    int blockSize = 0;
    double blockWeight = weightedLabels[0].second;
    for (unsigned int i = 0; i < weightedLabels.size(); i++) {
      if (std::abs(weightedLabels[i].second - blockWeight) < 1.0e-5) {
        blockSize++;
      }
      else {
        sizesSymmetryBlocks.push_back(blockSize);
        blockWeight = weightedLabels[i].second;
        blockSize = 1;
      }
    }
  }

  void NRG::reorderSymmetryBlocks() {
    int sum = 0;
    int passedCuts = 0;
    int nextCut = nrgInfo.statesKept + nrgInfo.statesAdded;
    for (int i = 0; i < (*this).getMaxBasisSize(); i++) {
      if (sum < nextCut && sum + sizesSymmetryBlocks[i] > nextCut) {
        int dummySum = sum;
        unsigned int j = i;
        while (dummySum != nextCut && j != sizesSymmetryBlocks.size()) {
          j++;
          if (dummySum + sizesSymmetryBlocks[j] <= nextCut) {
            dummySum += sizesSymmetryBlocks[j];
            moveBlock(j, i);
          }
        }
        if (j == sizesSymmetryBlocks.size()) {
          std::cerr << "Error whilst reordering the symmetry blocks to prevent ";
          std::cerr << "them from being introduced in different RG-steps. ";
          std::cerr << "The current routine cannot find a suitable block at ";
          std::cerr << "for one of the cuts, create more states in the scanning ";
          std::cerr << "routine or consider the block sizes and consider their ";
          std::cerr << "compatibily to the truncation points. " << std::endl;
          exit (-1);
        }
      }
      sum += sizesSymmetryBlocks[i];
      if (sum == nextCut) {
        passedCuts++;
        nextCut += nrgInfo.statesAdded;
      }
    }
  }

  void NRG::moveBlock (int originalPosition, int finalPosition) {
    int startBlock = std::accumulate(sizesSymmetryBlocks.begin(),
        sizesSymmetryBlocks.begin() + originalPosition, 0);
    std::vector<std::pair<std::string, double> > copiedBlock;
    for (int i = startBlock; i < startBlock +
        sizesSymmetryBlocks[originalPosition]; i++) {
      copiedBlock.push_back(weightedLabels[i]);
    }
    weightedLabels.erase(weightedLabels.begin() + startBlock,
        weightedLabels.begin() + startBlock +
        sizesSymmetryBlocks[originalPosition]);
    int indexInsertion = std::accumulate(sizesSymmetryBlocks.begin(),
        sizesSymmetryBlocks.begin() + finalPosition, 0);
    for (unsigned int i = 0; i < copiedBlock.size(); i++) {
      weightedLabels.insert(weightedLabels.begin() + indexInsertion,
          copiedBlock[i]); 
    }
    int temporary = sizesSymmetryBlocks[originalPosition];
    sizesSymmetryBlocks[originalPosition] = sizesSymmetryBlocks[finalPosition];
    sizesSymmetryBlocks[finalPosition] = temporary;
  }

}


