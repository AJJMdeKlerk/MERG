#ifndef NRG_H
#define NRG_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {

  class NRG {
    private:
      // Contains the required info for the routine
      QuenchInfo quenchInfo;
      ScanInfo scanInfo;
      NRGInfo nrgInfo;
      TimeEvolutionInfo timeEvolutionInfo;

      // Computed during the routines
      int currentNRGStep;
      int currentBasisSize;
      std::vector<std::pair<std::string, double> > weightedLabels;
      std::vector<int> sizesSymmetryBlocks;
      std::vector<LiebLinState> computationalBasis;
      Eigen::VectorXd energiesComputationalBasis;
      Eigen::VectorXd approxEigenstateEnergies;
      Eigen::MatrixXd approxEigenstateExpansions;
      Eigen::VectorXi selectedEigenstates;
      Eigen::MatrixXd selectedApproxEigenstateExpansions;
      Eigen::MatrixXd perturbingOperator;
      Eigen::MatrixXd nrgHamiltonian;

      // Output files
      std::string filenameEnergyOutput;
      std::string filenameReturnAmplitudeAndFidelityOutput;
      std::string filenameTimeEvolutionSelectedOperator;
      void setFilenameEnergyOutput();
      void setFilenameReturnAmplitudeAndFidelityOutput();
      void setFilenameTimeEvolutionSelectedOperator();

      // General preparation for all NRG routines
      bool hasFinishedBefore();
      void loadComputationalBasis();
      void checkSizeWeightedLabels();

      // Initialising the necessary objects and operators
      void initEnergiesComputationalBasis();
      void initApproxEigenstateEnergies();
      void initApproxEigenstateExpansions();
      void initSelectedApproxEigenstates();
      void initNRGHamiltonian();
      void initPerturbingOperator();

      // Updating the objects and operators at the start of a new NRG-step
      void updatePerturbingOperator();
      void updateNRGHamiltonian();
      void updateSelectedApproxEigenstates();
      void addNewStatesToSelectedApproxEigenstates();
      void computeSelectedApproxEigenstateExpansions();
      void computeWeightsAndSelectApproxEigenstates();

      // Processing the NRG-Hamiltonian of a given NRG-step
      void processNRGHamiltonian();
      void diagonaliseHamiltonianAndUpdateApproxEigenvectors();
      void writeResultsDiagonalisationToTerminal();
      void writeResultsDiagonalisationToFile();

      // Time evolution
      void computeReturnAmplitudeAndFidelity();
      void computeApproximateTimeEvolutionSelectedOperator();
      Eigen::MatrixXd getWeightedMatrixElementsTimeEvolution() const;
      std::vector<std::vector<int> >
        getIndicesSelectedSummandsTimeEvolution(const Eigen::MatrixXd&
            weightedMatrixElementsTimeEvolution) const;

      // Auxillliary
      double getMEPerturbingOperator(const int& i, const int& j) const;
      void computeSizesSymmetryBlocks();
      void reorderSymmetryBlocks();
      void moveBlock (int originalPosition, int finalPosition);

    public:
      std::string getFilenameEnergyOutput() const;
      std::string getFilenameReturnAmplitudeAndFidelityOutput() const;
      std::string getFilenameTimeEvolutionSelectedOperator() const;
      double getTargetEnergy() const;
      int getMaxNRGSteps() const;
      int getMaxBasisSize() const;
      int getCurrentNRGStep() const;
      int getCurrentBasisSize() const;

      NRG(const QuenchInfo& inputQuenchInfo, const ScanInfo& inputScanInfo,
          const NRGInfo& inputNRGInfo, const TimeEvolutionInfo&
          inputTimeEvolutionInfo);

  };

}

#endif
