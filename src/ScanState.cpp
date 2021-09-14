#include "ABACUS_VB.h"

namespace ABACUS_VB {

  std::queue<ScanState> ScanState::generateEqualMomentumDescendents() {
    std::queue<ScanState> descendents;
    std::queue<ScanState> auxilliaryDescendents;
    (*this).generateLeftMovers(auxilliaryDescendents);
    while(!auxilliaryDescendents.empty()) {
      auxilliaryDescendents.front().generateRightMovers(descendents);
      auxilliaryDescendents.pop();
    }
    return descendents;
  }

  std::queue<ScanState> ScanState::generateDescendents() {
    std::queue<ScanState> descendents;
    (*this).generateRightMovers(descendents);
    if (!rightMoversPresent)
      (*this).generateLeftMovers(descendents);
    return descendents;
  }

  void ScanState::generateLeftMovers(std::queue<ScanState>& descendents) {
    for (int i = indexRightmostLeftmover; i < doubledQuantumNumbs.size(); i++) {
      ScanState descendent(*this);
      if(descendent.leftmoveDoubledQuantumNumber(i))
        descendents.push(descendent);
    }
  }

  void ScanState::generateRightMovers(std::queue<ScanState>& descendents) {
    for (int i = 0; i <= indexLeftmostRightmover; i++) {
      ScanState descendent(*this);
      if(descendent.rightmoveDoubledQuantumNumber(i))
        descendents.push(descendent);
    }
  }

  bool ScanState::leftmoveDoubledQuantumNumber(const int& i) {
    if ((moverTypes(i) != 1) && (i == 0
          || doubledQuantumNumbs(i) - 2 != doubledQuantumNumbs(i - 1))) {
      doubledQuantumNumbs(i) -= 2;
      moverTypes(i) = -1;
      indexRightmostLeftmover = i;
      return true;
    }
    return false;
  }

  bool ScanState::rightmoveDoubledQuantumNumber(const int& i) {
    if ((moverTypes(i) != -1) && (i == doubledQuantumNumbs.size() - 1
          || doubledQuantumNumbs(i) + 2 != doubledQuantumNumbs(i + 1))) {
      doubledQuantumNumbs(i) += 2;
      moverTypes(i) = 1;
      rightMoversPresent = true;
      indexLeftmostRightmover = i;
      return true;
    }
    return false;
  }

  void ScanState::findIndexLeftmostRightmover() {
    bool foundLeftmostRightmover = false;
    int i = 0;
    while(!foundLeftmostRightmover && i != doubledQuantumNumbs.size()) {
      if (moverTypes(i) == 1) {
        indexLeftmostRightmover = i;
        foundLeftmostRightmover = true;
      }
      i++;
    }
    if (!foundLeftmostRightmover)
      indexLeftmostRightmover = doubledQuantumNumbs.size() - 1;
  }

  void ScanState::findIndexRightmostLeftmover() {
    bool foundRightmostLeftmover = false;
    int i = doubledQuantumNumbs.size() - 1;
    while(!foundRightmostLeftmover && i != -1) {
      if (moverTypes(i) == -1) {
        indexRightmostLeftmover = i;
        foundRightmostLeftmover = true;
      }
      i--;
    }
    if (!foundRightmostLeftmover)
      indexRightmostLeftmover = 0;
  }

  void ScanState::findMoverTypes() {
    moverTypes = Eigen::VectorXi::Zero(doubledQuantumNumbs.size());
    for (int i = 0; i < doubledQuantumNumbs.size(); i++) {
      if (doubledQuantumNumbs(i) < doubledQuantumNumbsSeedState(i))
        moverTypes(i) = -1;
      else if (doubledQuantumNumbs(i) > doubledQuantumNumbsSeedState(i))
        moverTypes(i) = 1;
    }
    foundMoverTypes = true;
  }

  double getWeight(const std::string& labelAndWeight) {
    std::stringstream splittableString(labelAndWeight);
    std::string label;
    splittableString >> label;
    double weight;
    splittableString >> weight;
    return weight;
  }

  std::string getLabel(const std::string& labelAndWeight) {
    std::stringstream splittableString(labelAndWeight);
    std::string label;
    splittableString >> label;
    return label;
  }

  bool isApprovedByWeight(const ScanInfo& scanInfo, const double& weight) {
    if (std::abs(weight) <= scanInfo.maxAbsWeight)
      return true;
    else
      return false;
  }

  bool isApprovedByWeight(const ScanInfo& scanInfo, const std::string& output) {
    std::stringstream readableOutput(output);
    std::string label;
    double weight;
    readableOutput >> label;
    readableOutput >> weight;
    if (std::abs(weight) <= scanInfo.maxAbsWeight)
      return true;
    else
      return false;
  }

  bool ScanState::isApprovedByWeightAndMomentum(const ScanInfo& scanInfo) {
    double weight = (*this).getWeight(scanInfo);
    int differenceSumQuantumNumbs = (*this).getSumDoubledQuantumNumbs() -
      scanInfo.seedState.getSumDoubledQuantumNumbs();
    if (differenceSumQuantumNumbs >= scanInfo.changeDoubledQuantumNumbs.first &&
        differenceSumQuantumNumbs <= scanInfo.changeDoubledQuantumNumbs.second &&
        weight <= scanInfo.maxAbsWeight)
      return true;
    else
      return false;
  }

  bool ScanState::isApprovedByMomentum(const ScanInfo& scanInfo) {
    int differenceSumQuantumNumbs = (*this).getSumDoubledQuantumNumbs() -
      scanInfo.seedState.getSumDoubledQuantumNumbs();
    if (differenceSumQuantumNumbs >= scanInfo.changeDoubledQuantumNumbs.first &&
        differenceSumQuantumNumbs <= scanInfo.changeDoubledQuantumNumbs.second)
      return true;
    else
      return false;
  }

  void ScanState::resetMoverInfo() {
    moverTypes.setZero();
    foundMoverTypes = true;
    int numberOfParticles = doubledQuantumNumbs.size();
    indexLeftmostRightmover = numberOfParticles - 1;
    indexRightmostLeftmover = 0;
    rightMoversPresent = false;
  }

  Eigen::VectorXi convertAbsoluteLabelToDoubledQuantumNumbs
    (const std::string& absoluteLabel, const ScanInfo& scanInfo) {
      Eigen::VectorXi doubledQuantumNumbs =
        scanInfo.seedState.getDoubledQuantumNumbs();

    std::stringstream splittableLabel(absoluteLabel);
    int firstDoubledQuantumNumber;
    std::string label;
    char deliminator;
    splittableLabel >> firstDoubledQuantumNumber;
    splittableLabel >> deliminator;
    splittableLabel >> label;
    std::string binaryLabel = convertHexToBinary(label);
    int distanceBetweenQuantumNumbs = 0;
    int quantumNumbsSet = 0;
    int i = 0;
    while(quantumNumbsSet < scanInfo.seedState.getParticleNumber()) {
      if(binaryLabel.at(i) == '0')
        distanceBetweenQuantumNumbs++;
      else if (binaryLabel.at(i) == '1' && quantumNumbsSet == 0) {
        doubledQuantumNumbs(0) = firstDoubledQuantumNumber;
        distanceBetweenQuantumNumbs = 0;
        quantumNumbsSet++;
      }
      else if (binaryLabel.at(i) == '1' && quantumNumbsSet != 0) {
        doubledQuantumNumbs(quantumNumbsSet) =
          doubledQuantumNumbs(quantumNumbsSet - 1) + 2 *
          (distanceBetweenQuantumNumbs + 1);
        distanceBetweenQuantumNumbs = 0;
        quantumNumbsSet++;
      }
      else {
        std::cerr << "Something went wrong in converting the binary label. ";
        std::cerr << "Exiting." << std::endl;
        exit(-1);
      }
     i++;
    }
    return doubledQuantumNumbs;
  }

  std::string convertBinaryToHex(const std::string& binaryNumber) {
    std::string paddedBinaryNumber = binaryNumber;
    while (paddedBinaryNumber.length() % 4 != 0) {
      paddedBinaryNumber = "0" + paddedBinaryNumber;
    }
    std::stringstream hexNumber;
    std::stringstream hexDigit;
    int i = 0;
    while (hexNumber.str().length() != paddedBinaryNumber.length() / 4) {
      hexDigit.str("");
      hexDigit << paddedBinaryNumber.at(i);
      hexDigit << paddedBinaryNumber.at(i+1);
      hexDigit << paddedBinaryNumber.at(i+2);
      hexDigit << paddedBinaryNumber.at(i+3);
      if (hexDigit.str() == "0000")
        hexNumber << '0';
      else if (hexDigit.str() == "0001")
        hexNumber << '1';
      else if (hexDigit.str() == "0010")
        hexNumber << '2';
      else if (hexDigit.str() == "0011")
        hexNumber << '3';
      else if (hexDigit.str() == "0100")
        hexNumber << '4';
      else if (hexDigit.str() == "0101")
        hexNumber << '5';
      else if (hexDigit.str() == "0110")
        hexNumber << '6';
      else if (hexDigit.str() == "0111")
        hexNumber << '7';
      else if (hexDigit.str() == "1000")
        hexNumber << '8';
      else if (hexDigit.str() == "1001")
        hexNumber << '9';
      else if (hexDigit.str() == "1010")
        hexNumber << 'a';
      else if (hexDigit.str() == "1011")
        hexNumber << 'b';
      else if (hexDigit.str() == "1100")
        hexNumber << 'c';
      else if (hexDigit.str() == "1101")
        hexNumber << 'd';
      else if (hexDigit.str() == "1110")
        hexNumber << 'e';
      else if (hexDigit.str() == "1111")
        hexNumber << 'f';
      i += 4;
    }
    return hexNumber.str();
  }

  std::string convertHexToBinary(const std::string& hexNumber) {
    std::stringstream binaryNumber;
    char digitHexNumber;
    for (unsigned int i = 0; i < hexNumber.length(); i++) {
      digitHexNumber = hexNumber.at(i);
      if (digitHexNumber == '0')
        binaryNumber << "0000";
      else if (digitHexNumber == '1')
        binaryNumber << "0001";
      else if (digitHexNumber == '2')
        binaryNumber << "0010";
      else if (digitHexNumber == '3')
        binaryNumber << "0011";
      else if (digitHexNumber == '4')
        binaryNumber << "0100";
      else if (digitHexNumber == '5')
        binaryNumber << "0101";
      else if (digitHexNumber == '6')
        binaryNumber << "0110";
      else if (digitHexNumber == '7')
        binaryNumber << "0111";
      else if (digitHexNumber == '8')
        binaryNumber << "1000";
      else if (digitHexNumber == '9')
        binaryNumber << "1001";
      else if (digitHexNumber == 'a')
        binaryNumber << "1010";
      else if (digitHexNumber == 'b')
        binaryNumber << "1011";
      else if (digitHexNumber == 'c')
        binaryNumber << "1100";
      else if (digitHexNumber == 'd')
        binaryNumber << "1101";
      else if (digitHexNumber == 'e')
        binaryNumber << "1110";
      else if (digitHexNumber == 'f')
        binaryNumber << "1111";
    }
    return binaryNumber.str();
  }

  Eigen::VectorXi convertLabelToDoubledQuantumNumbs(const std::string& label,
      const ScanInfo& scanInfo) {
    Eigen::VectorXi doubledQuantumNumbsSeedState =
      scanInfo.seedState.getDoubledQuantumNumbs();
    Eigen::VectorXi doubledQuantumNumbs = doubledQuantumNumbsSeedState;
    std::stringstream splittableLabel(label);
    int doubledQuantumNumber;
    char deliminator;
    if (splittableLabel.peek() == '|')
      splittableLabel >> deliminator;
    while (splittableLabel >> doubledQuantumNumber) {
      auto i = 0;
      while (i < doubledQuantumNumbs.size()) {
        if (doubledQuantumNumber == doubledQuantumNumbsSeedState(i))
          break;
        i++;
      }
      splittableLabel >> deliminator;
      splittableLabel >> doubledQuantumNumber;
      doubledQuantumNumbs(i) = doubledQuantumNumber;
      splittableLabel >> deliminator;
    }
    return doubledQuantumNumbs;
  }

}
