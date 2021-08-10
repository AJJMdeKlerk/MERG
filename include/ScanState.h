#ifndef ScanState_H
#define ScanState_H

#include "ABACUS_VB.h"

namespace ABACUS_VB {

  class LiebLinState;

  class ScanState {
    private:
      Eigen::VectorXi doubledQuantumNumbsSeedState;
      Eigen::VectorXi doubledQuantumNumbs;
      Eigen::VectorXi moverTypes;

      bool foundMoverTypes;
      bool rightMoversPresent;
      int indexLeftmostRightmover;
      int indexRightmostLeftmover;

      void findIndexLeftmostRightmover();
      void findIndexRightmostLeftmover();
      void findMoverTypes();

      void generateLeftMovers(std::queue<ScanState>& descendents);
      void generateRightMovers(std::queue<ScanState>& descendents);
    public:
      std::queue<ScanState> generateDescendents();
      std::queue<ScanState> generateEqualMomentumDescendents();
      bool leftmoveDoubledQuantumNumber(const int& i);
      bool rightmoveDoubledQuantumNumber(const int& i);
      void resetMoverInfo();

      std::string getLabelAndWeight(const ScanInfo& scanInfo) const;
      std::string getAdditionalOutput(const ScanInfo& scanInfo) const;
      std::string getLabel() const;
      std::string getAbsoluteLabel() const;
      double getWeight(const ScanInfo& scanInfo) const;
      int getSumDoubledQuantumNumbs() const;
      int getNumberOfDoubledQuantumNumbs() const;
      Eigen::VectorXi getDoubledQuantumNumbs() const;
      Eigen::VectorXi getDoubledQuantumNumbsSeedState() const;
      Eigen::VectorXi getMoverTypes() const;

      void setFromScanInfo(const ScanInfo& scanInfo);
      bool isApprovedByWeightAndMomentum(const ScanInfo& scanInfo);
      bool isApprovedByMomentum(const ScanInfo& scanInfo);

      ScanState();
      ScanState(const ScanInfo& scanInfo);
      ScanState(const ScanInfo& scanInfo, const LiebLinState& state);
      ScanState(const ScanInfo& scanInfo, const std::string& label);
      ScanState(const ScanInfo& scanInfo,
          const Eigen::VectorXi& inputDoubledQuantumNumbs);
  };

  bool isApprovedByWeight(const ScanInfo& scanInfo, const std::string&
      output);
  bool isApprovedByWeight(const ScanInfo& scanInfo, const double& weight);
  double getWeight(const std::string& labelAndWeight);
  std::string getLabel(const std::string& labelAndWeight);
  Eigen::VectorXi convertLabelToDoubledQuantumNumbs(const std::string& label,
      const ScanInfo& scanInfo);
  Eigen::VectorXi convertAbsoluteLabelToDoubledQuantumNumbs(const
      std::string& binaryBasedLabel, const ScanInfo& scanInfo);
  std::string convertHexToBinary(const std::string& hexNumber);
  std::string convertBinaryToHex(const std::string& hexNumber);


}

#endif
