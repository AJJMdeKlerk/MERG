#include "ABACUS_VB.h"

namespace ABACUS_VB {
 
  void LiebLinState::setDoubledQuantumNumbs(const Eigen::VectorXi& input) {
    doubledQuantumNumbs = input;
    (*this).computeAll();
  }

  std::string LiebLinState::getInfoForFilename() const {
    std::stringstream info;
    mkdir("data/", S_IRUSR | S_IWUSR | S_IXUSR);
    info << "data/";
    info << "N" << (*this).getParticleNumber() << "_";
    info << "L" << (*this).getSystemLength() << "_";
    info << "C" << (*this).getInteractionStrength() << "/";
    mkdir(info.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    info << "seedState" << doubledQuantumNumbs.transpose().format(CSVReduced);
    info << "/";
    mkdir(info.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    return info.str();
  }

  double LiebLinState::getRapidity(int i) const {
    if (!allComputed) {
      std::cerr << "A rapidity was requested before being computed. ";
      std::cerr << "Exiting program." << std::endl;
      exit(-1);
    }
    return rapidities(i);
  }

  Eigen::VectorXd LiebLinState::getRapidities() const {
    if (!allComputed) {
      std::cerr << "The rapidities were requested before being computed. ";
      std::cerr << "Exiting program." << std::endl;
      exit(-1);
    }
    return rapidities.transpose();
  }

  Eigen::MatrixXd LiebLinState::getGaudinMatrix() const {
    if (!allComputed) {
      std::cerr << "The Gaudin matrix were requested before being computed. ";
      std::cerr << "Exiting program." << std::endl;
      exit(-1);
    }
    return gaudinMatrix;
  }

  double LiebLinState::getEnergy() const {
    if (!allComputed) {
      std::cerr << "The energy was requested before being computed.";
      std::cerr << "Exiting program." << std::endl;
      exit(-1);
    }
    return energy;
  }

  double LiebLinState::getMomentum() const {
    if (!allComputed) {
      std::cerr << "The momentum was requested before being computed.";
      std::cerr << "Exiting program." << std::endl;
      exit(-1);
    }
    return momentum;
  }

  double LiebLinState::getThirdMoment() const {
    if (!allComputed) {
      std::cerr << "The third moment was requested its computation.";
      std::cerr << "Exiting program." << std::endl;
      exit(-1);
    }
    return thirdMoment;
  }

  double LiebLinState::getLnNorm() const {
    if (!allComputed) {
      std::cerr << "The norm was requested before being computed.";
      std::cerr << "Exiting program." << std::endl;
      exit(-1);
    }
    return lnNorm;
  }

  SystemInfo LiebLinState::getSystemInfo() const {
    return systemInfo;
  }

  Eigen::VectorXi LiebLinState::getDoubledQuantumNumbs() const {
    return doubledQuantumNumbs;
  }

  int LiebLinState::getSumDoubledQuantumNumbs() const {
    return doubledQuantumNumbs.sum();
  }

  SystemInfo LiebLinState::getsystemInfo() const {
    return systemInfo;
  }

  int LiebLinState::getParticleNumber() const {
    return systemInfo.particleNumber;
  }

  double LiebLinState::getInteractionStrength() const {
    return systemInfo.interactionStrength;
  }

  double LiebLinState::getSystemLength() const {
    return systemInfo.systemLength;
  }

  double LiebLinState::getMaxErrorRapidities() const {
    return maxErrorRapidities;
  }

  int LiebLinState::getMaxNewtonIterations() const {
    return maxNewtonIterations;
  }

}
