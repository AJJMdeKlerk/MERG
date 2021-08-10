#include "ABACUS_VB.h"

namespace ABACUS_VB {
  
  void LiebLinState::computeAll() {
    (*this).computeRapidities();
    (*this).computeLnNorm();
    (*this).computeMoments();
    allComputed = true;
  }

  void LiebLinState::computeRapidities() {
    // Solves the Bethe equations up to a given accuracy using the
    // Newton-Rhapson algorithm. The Jacobian is the Gaudin matrix.
    Eigen::VectorXd changeTrialSolution((*this).getParticleNumber());
    Eigen::VectorXd deviationRapidities((*this).getParticleNumber());
    deviationRapidities = -(*this).BAEs(rapidities);
    int newtonIter = 0;
    while (deviationRapidities.norm() >= maxErrorRapidities 
      && newtonIter < maxNewtonIterations) {
      (*this).computeGaudinMatrix();
      deviationRapidities = -(*this).BAEs(rapidities);
      changeTrialSolution = (gaudinMatrix.partialPivLu()).solve
        (deviationRapidities);
      rapidities += changeTrialSolution; 
      newtonIter++;
    }
    if (newtonIter == maxNewtonIterations) {
      std::cout << "The Newton Rhapson algorithm used to solve the ";
      std::cout << "logarithmic Bethe equations did not converge. ";
      std::cout << "Exciting to prevent the possible generation of ";
      std::cout << "inaccurate data." << std::endl;
      exit(0);
    }
  }

  void LiebLinState::computeGaudinMatrix() {
    gaudinMatrix = Eigen::MatrixXd::Zero
      ((*this).getParticleNumber(), (*this).getParticleNumber());
    tbb::parallel_for(0, (*this).getParticleNumber(),
      [&](int j) {
        for (int i = 0; i < (*this).getParticleNumber(); i++){
          if (i == j) {
            gaudinMatrix(i, j) += (*this).getSystemLength() - (*this).kernel(i, j);
            for (int k = 0; k < (*this).getParticleNumber(); k++) {
              gaudinMatrix(i, j) += (*this).kernel(i, k);
            }
          }
          else {
            gaudinMatrix(i, j) -= (*this).kernel(i, j);
          }
        }
      }
    );
  }

  void LiebLinState::computeLnNorm() {
    // The norm of a Bethe state is equal to the determinant of the 
    // Gaudin matrix, but here we include an additional part in the 
    // norm as is customary.
  	(*this).computeGaudinMatrix();
    lnNorm = std::real(lnDetInvertibleMatrix(gaudinMatrix));
    lnNorm += (*this).getParticleNumber()
      * std::log((*this).getInteractionStrength());
  	for (int i = 0; i < (*this).getParticleNumber(); i++) {
  		for (int j = i + 1; j < (*this).getParticleNumber(); j++) {
  			lnNorm += std::log(1 + std::pow((*this).getInteractionStrength(), 2.0) 
            / std::pow(rapidities(i) - rapidities(j), 2.0));
  		}
  	}
  }

  void LiebLinState::computeMoments() {
    // The moments we are interested in, such as the energy, momentum and third
    // moment necessary for the g2 matrix elements, can be obtained from the
    // rapidities
    momentum = 0.0;
    energy = 0.0;
    thirdMoment = 0.0;
    for (int i = 0; i < (*this).getParticleNumber(); i++) {
      momentum += rapidities(i);
      energy += std::pow(rapidities(i), 2.0);
      thirdMoment += std::pow(rapidities(i), 3.0);
    }
  }

  Eigen::VectorXd LiebLinState::BAEs(const Eigen::VectorXd& trialSolution) {
    // Gives 0 if the input satisfies the (logarithmic) Bethe equations
    if (trialSolution.size() != (*this).getParticleNumber()) {
      std::cerr << "Error, the size of the trial solution plugged into the ";
      std::cerr << "Bethe equations is not the same size as the number of ";
      std::cerr << "particles.";
      exit(-1);
    }
    Eigen::VectorXd result = Eigen::VectorXd::Zero
      ((*this).getParticleNumber());
    tbb::parallel_for(0, (*this).getParticleNumber(),
      [&](int i){
        for (int j = 0; j < (*this).getParticleNumber(); j++) {
          result(i) += 2.0 * std::atan((trialSolution(i) 
            - trialSolution(j)) / (*this).getInteractionStrength());
        }
        result(i) += trialSolution(i) * (*this).getSystemLength();
        result(i) -= M_PI * doubledQuantumNumbs(i);
      }
    );
    return result;
  }

  Eigen::VectorXd LiebLinState::interactionDerivativeBAEs
    (const Eigen::VectorXd& trialSolution) {
    // Gives 0 if the input satisfies the derivative with respect to the 
    // interaction strength of the (logarithmic) Bethe equations
    if (trialSolution.size() != (*this).getParticleNumber()) {
      std::cerr << "Error, the size of the trial solution plugged into the ";
      std::cerr << "derivative of the Bethe equations w.r.t. the interaction ";
      std::cerr << "strength is not the same size as the number of particles.";
      exit(-1);
    }
    Eigen::VectorXd result = Eigen::VectorXd::Zero
      ((*this).getParticleNumber());
    tbb::parallel_for(0, (*this).getParticleNumber(),
      [&](int i){
        for (int j = 0; j < (*this).getParticleNumber(); j++) {
          result(i) -= (1.0 / (*this).getInteractionStrength())
            * (rapidities(i) - rapidities(j))
            * (*this).kernel(i, j);
          result(i) += (*this).kernel(i, j) 
            * (trialSolution(i) - trialSolution(j));
        }
        result(i) += (*this).getSystemLength() * trialSolution(i);
      }
    );
    return result;
  }

  double LiebLinState::kernel(const int i, const int j) const {
  	return (2.0 * (*this).getInteractionStrength()) 
      / (std::pow(rapidities(i) - rapidities(j), 2.0)
        + std::pow((*this).getInteractionStrength(), 2.0));
  }

  void LiebLinState::displayState() const {
    if (!allComputed) {
      std::cerr << "The state was requested to be displayed before it was ";
      std::cerr << "computed. Exiting program." << std::endl;
      exit(-1);
    }
    else {
      std::cout << "Doubled quantum numbers: ";
      std::cout << (*this).getDoubledQuantumNumbs().transpose() << std::endl;
      std::cout << "Rapidities: ";
      std::cout << (*this).getRapidities().transpose() << std::endl;
      std::cout << "Energy: ";
      std::cout << (*this).getEnergy() << std::endl;
      std::cout << std::endl;
    }
  }

  void LiebLinState::checkPositiveInteractionStrength() {
    if ((*this).getInteractionStrength() <= 0.0) {
  		std::cout << "A non-positive value for the interaction parameter ";
      std::cout << "is not supported in this version. " << std::endl;
  		exit(-1);
  	}
  }

}
