#include "ABACUS_VB.h"

namespace ABACUS_VB {
  
  double rho(const LiebLinState& leftState, const LiebLinState& rightState) {
    if(!isCompatible(leftState, rightState)){
      std::cerr << "Error, you attempted to compute the matrix element of ";
      std::cerr << "the density operator for two states which are not ";
      std::cerr << "compatible. " << std::endl;
      exit(-1);
    }
    if(leftState.getDoubledQuantumNumbs() ==
        rightState.getDoubledQuantumNumbs()) {
      return leftState.getParticleNumber() / leftState.getSystemLength();
    }
    else {
      return std::real(std::exp(lnRhoOffDiagonal(leftState, rightState)));
    }
  }

  std::complex<double> lnRhoOffDiagonal(const LiebLinState& leftState, const
      LiebLinState& rightState) { 

    double J = rightState.getMomentum() - leftState.getMomentum();

    std::complex<double> result = std::log(std::complex<double>(II * J)); 
    for(int i = 0; i < leftState.getParticleNumber(); i++) {
      for(int j = 0; j < leftState.getParticleNumber(); j++) {
        result += std::log((rightState.getRapidity(i) -
            rightState.getRapidity(j) + II *
            rightState.getInteractionStrength())
          / (leftState.getRapidity(i) - rightState.getRapidity(j)));
      }
    }

    for(int i = 0; i < leftState.getParticleNumber(); i++) {
      result += std::log(2.0 * II * std::imag(std::exp(
              lnVPlus(leftState.getRapidities(), rightState.getRapidities(), 
          leftState.getInteractionStrength(), i))));
    }

    // Arbitrary number in the expression
    double lambdaP = M_PI;
    
    // Construct the matrix whose determinant will be added to the result
    // We start with the matrix U
    Eigen::MatrixXd uMatrix = Eigen::MatrixXd::Ones 
      (leftState.getParticleNumber(), leftState.getParticleNumber());
    for (int i = 0; i < leftState.getParticleNumber(); i++) {
      std::complex<double> lnProduct = 0.0;
      for (int m = 0; m < leftState.getParticleNumber(); m++) {
        if (m != i) {
          lnProduct += std::log(std::complex<double>((leftState.getRapidity(m) 
            - rightState.getRapidity(i)) / (rightState.getRapidity(m) 
            - rightState.getRapidity(i))));
        }
        else {
          lnProduct += std::log(std::complex<double>(leftState.getRapidity(m) 
            - rightState.getRapidity(i)));
        }
      }
      uMatrix.row(i) *= std::real(std::exp(lnProduct)) / (2.0 *
        std::imag(std::exp(lnVPlus(leftState.getRapidities(),
          rightState.getRapidities(), leftState.getInteractionStrength(), i))));
      for (int j = 0; j < leftState.getParticleNumber(); j++) {
        uMatrix(i, j) *= (rightState.kernel(i, j) -
            kernel(lambdaP, rightState.getRapidity(j),
              rightState.getInteractionStrength()));
      }
    }

    Eigen::MatrixXd fullMatrix = Eigen::MatrixXd::Identity
      (leftState.getParticleNumber(), leftState.getParticleNumber());
    fullMatrix += uMatrix;
    result += lnDetInvertibleMatrix(fullMatrix);

    // Divide by term with arbirary V's
    result -= std::log(std::complex<double>(2.0 * II * std::imag(std::exp(
              lnVPlusArbitrary(leftState.getRapidities(), 
                rightState.getRapidities(), lambdaP, 
                rightState.getInteractionStrength())))));

    // Divide by the norms
    result -= 0.5 * leftState.getLnNorm();
    result -= 0.5 * rightState.getLnNorm();

    if (std::isnan(std::real(std::exp(result)))) {
      std::cerr << "Error, computing an off-diagonal g2 matrix element ";
      std::cerr << "returned a NaN. Exiting." << std::endl;
      std::cerr << "The states it failed had Ix2s given by ";
      std::cerr << leftState.getDoubledQuantumNumbs().transpose() << std::endl;
      std::cerr << rightState.getDoubledQuantumNumbs().transpose() << std::endl;
      exit(-1);
    }
    else 
      return result;
  }

  double G2(const LiebLinState& leftState, const LiebLinState& rightState) {
    if(!isCompatible(leftState, rightState)){
      std::cerr << "Error, you attempted to compute the g2 matrix element of ";
      std::cerr << "two states which are not compatible. " << std::endl;
      exit(-1);
    }
    if(leftState.getDoubledQuantumNumbs() == 
        rightState.getDoubledQuantumNumbs()){
      return G2Diagonal(leftState);
    }
    else {
      return std::real(std::exp(lnG2OffDiagonal(leftState, rightState)));
    }
  }

  double G2Diagonal(const LiebLinState& state) {
    double result = 0.0;
    Eigen::VectorXd interactionDerivativeRapidities = 
      computeInteractionDerivativeRapidities(state);
    result = 2.0 * state.getRapidities().dot(interactionDerivativeRapidities);
    result /= state.getSystemLength();
    return result;
  }

  Eigen::VectorXd computeInteractionDerivativeRapidities 
    (const LiebLinState& state) {
    Eigen::VectorXd trialSolution = initInteractionDerivativeRapidities(state);
    Eigen::VectorXd deviationTrialSolution = 
      - interactionDerivativeBAEs(state, trialSolution);
    Eigen::VectorXd changeTrialSolution;  
    int newtonIter = 0;
    while (deviationTrialSolution.norm() >= state.getMaxErrorRapidities()
      && newtonIter < state.getMaxNewtonIterations()) {
      deviationTrialSolution = 
        - interactionDerivativeBAEs(state, trialSolution);
      changeTrialSolution = (state.getGaudinMatrix().partialPivLu()).solve
        (deviationTrialSolution);
      trialSolution += changeTrialSolution; 
      newtonIter++;
    }
    if (newtonIter == state.getMaxNewtonIterations() - 1) {
      std::cerr << "The Newton Rhapson algorithm used to solve the ";
      std::cerr << "derivative of the logarithmic Bethe equations with";
      std::cerr << "respect to the interaction strength did not converge. ";
      std::cerr << "Exciting to prevent the possible generation of ";
      std::cerr << "inaccurate data." << std::endl;
      exit(0);
    }
    return trialSolution;
  }

  Eigen::VectorXd interactionDerivativeBAEs 
    (const LiebLinState& state, const Eigen::VectorXd& trialSolution) {
    // Gives 0 if the input satisfies the derivative with respect to the 
    // interaction strength of the (logarithmic) Bethe equations
    Eigen::VectorXd result = Eigen::VectorXd::Zero(state.getParticleNumber());
    tbb::parallel_for(0, state.getParticleNumber(),
      [&](int i){
        for (int j = 0; j < state.getParticleNumber(); j++) {
          result(i) -= (1.0 / state.getInteractionStrength())
            * (state.getRapidity(i) - state.getRapidity(j))
            * state.kernel(i, j);
          result(i) += state.kernel(i, j) 
            * (trialSolution(i) - trialSolution(j));
        }
        result(i) += state.getSystemLength() * trialSolution(i);
      }
    );
    return result;
  }

  Eigen::VectorXd initInteractionDerivativeRapidities 
    (const LiebLinState& state) {
    Eigen::VectorXd result(state.getParticleNumber());
    for (int i = 0; i < state.getParticleNumber(); i++) {
      for (int k = 0; k < state.getParticleNumber(); k++) {
        result(i) = (1.0 / state.getInteractionStrength())
          * state.kernel(i, k) * state.kernel(i, k);
      }
    }
    return result;
  }

  std::complex<double> lnG2OffDiagonal(const LiebLinState& leftState, 
      const LiebLinState& rightState) {
    // The expression implemented here is Eq. (22) of Piroli, L., & Calabrese,
    // P. (2015). Exact formulas for the form factors of local operators in the
    // Lieb–Liniger model. Journal of Physics A: Mathematical and Theoretical,
    // 48(45), 454002. 
    double J = 0.0; 
    J += std::pow(rightState.getMomentum() - leftState.getMomentum(), 4.0);
    J -= 4.0 * (rightState.getMomentum() - leftState.getMomentum()) 
      * (rightState.getThirdMoment() - leftState.getThirdMoment());
    J += 3.0 * std::pow(rightState.getEnergy() - leftState.getEnergy(), 2.0);

    // We add the prefactors not involving any products
    std::complex<double> result = std::log(std::complex<double>(J)); 
    result += std::log(std::complex<double> (std::pow(-1.0, 
      leftState.getParticleNumber())));
    result -= std::log(6.0 * leftState.getInteractionStrength());

    for(int i = 0; i < leftState.getParticleNumber(); i++) {
      for(int j = 0; j < leftState.getParticleNumber(); j++) {
        result += std::log(std::complex<double> ((rightState.getRapidity(i) -
            rightState.getRapidity(j) + II *
            rightState.getInteractionStrength())
          / (rightState.getRapidity(i) - leftState.getRapidity(j))));
      }
    }

    for(int i = 0; i < leftState.getParticleNumber(); i++) {
      result += std::log(2.0 * II * std::imag(std::exp(
              lnVPlus(leftState.getRapidities(), rightState.getRapidities(), 
          leftState.getInteractionStrength(), i))));
    }

    // Arbitrary numbers in the expression
    double lambdaP = rightState.getRapidity(0);
    double lambdaS = rightState.getRapidity(leftState.getParticleNumber() - 1);
    
    // Construct the matrix whose determinant will be added to the result
    // We start with the matrix U
    Eigen::MatrixXd uMatrix = Eigen::MatrixXd::Ones 
      (leftState.getParticleNumber(), leftState.getParticleNumber());
    for (int i = 0; i < leftState.getParticleNumber(); i++) {
      std::complex<double> lnProduct = 0.0;
      for (int m = 0; m < leftState.getParticleNumber(); m++) {
        if (m != i) {
          lnProduct += std::log(std::complex<double>((leftState.getRapidity(m) 
            - rightState.getRapidity(i)) / (rightState.getRapidity(m) 
            - rightState.getRapidity(i))));
        }
        else {
          lnProduct += std::log(std::complex<double>(leftState.getRapidity(m) 
            - rightState.getRapidity(i)));
        }
      }

      uMatrix.row(i) *= std::real(std::exp(lnProduct)) / (2.0 *
        std::imag(std::exp(lnVPlus(leftState.getRapidities(),
          rightState.getRapidities(), leftState.getInteractionStrength(), i))));

      for (int j = 0; j < leftState.getParticleNumber(); j++) {
         uMatrix(i,j) *= (kernel(rightState.getRapidity(i),
               rightState.getRapidity(j), rightState.getInteractionStrength())
                  - kernel(lambdaP, rightState.getRapidity(j),
                      rightState.getInteractionStrength())
                        * kernel(lambdaS, rightState.getRapidity(i),
                            rightState.getInteractionStrength()));
        if (i == j) {
          uMatrix(i,j) += 1.0;
        }
      }
    }
    result += lnDetInvertibleMatrix(uMatrix);

    // Divide by term with arbirary V's
    result -= std::log(std::complex<double>(2.0 * II * std::imag(std::exp(
              lnVPlusArbitrary(leftState.getRapidities(), 
                rightState.getRapidities(), lambdaP, 
                rightState.getInteractionStrength())))));
    result -= std::log(std::complex<double>(2.0 * II * std::imag(std::exp(
              lnVPlusArbitrary(leftState.getRapidities(), 
                rightState.getRapidities(), lambdaS, 
                rightState.getInteractionStrength())))));

    // Divide by the norms
    result -= 0.5 * leftState.getLnNorm();
    result -= 0.5 * rightState.getLnNorm();

    if (std::isnan(std::real(std::exp(result)))) {
      std::cerr << "Error, computing an off-diagonal g2 matrix element ";
      std::cerr << "returned a NaN. Exiting." << std::endl;
      std::cerr << "The states it failed had Ix2s given by: " << std::endl;
      std::cerr << leftState.getDoubledQuantumNumbs().transpose() << std::endl;
      std::cerr << rightState.getDoubledQuantumNumbs().transpose() << std::endl;
      exit(-1);
    }
    else 
      return result;
  }

  bool isCompatible(const LiebLinState& leftState, 
    const LiebLinState& rightState) {
    bool compatible = true;
    // Exits the program if we are trying to use non-compatible states
    if (leftState.getParticleNumber() != rightState.getParticleNumber()){
      std::cerr << "Error, the states used as input in lnG2 should have an ";
      std::cerr << "equal number of particles. Exciting. " << std::endl;
      compatible = false;
    }
    if (leftState.getInteractionStrength() 
        != rightState.getInteractionStrength()){
      std::cerr << "Error, the states used as input in lnG2 should have ";
      std::cerr << "the same interaction strength. Exciting. " << std::endl;
      compatible = false;
    }
    if (leftState.getSystemLength() != leftState.getSystemLength()){
      std::cerr << "Error, the states used as input in lnG2 should have ";
      std::cerr << "the same length. Exciting. " << std::endl;
      compatible = false;
    }
    return compatible;
  }

  double kernel(const double& lambda, const double& mu, const double& c) {
    return 2.0 * c / (std::pow(lambda - mu, 2.0) + std::pow(c, 2.0));
  }

  std::complex<double> lnVPlusArbitrary(const Eigen::VectorXd& mu, 
      const Eigen::VectorXd& lambda, std::complex<double> lambdaP, 
      double c) {
    if (mu.size() != lambda.size()) {
      std::cerr << "Error, the vectors serving as input to ";
      std::cerr << "should be of equal length. Exciting program. ";
      exit(-1);
    }
    std::complex<double> result = 0.0;
    for (int j = 0; j < mu.size(); j++) {
      result += std::log(std::complex<double> ((mu(j) - lambdaP + II * c)
          / (lambda(j) - lambdaP + II * c)));
    }
    return result;
  }

  std::complex<double> lnVPlus(const Eigen::VectorXd& mu, 
      const Eigen::VectorXd& lambda, double c, int i) {
    if (mu.size() != lambda.size()) {
      std::cerr << "Error, the vectors serving as input to ";
      std::cerr << "should be of equal length. Exciting program. ";
      exit(-1);
    }
    std::complex<double> result = 0.0;
    for (int j = 0; j < mu.size(); j++) {
      result += std::log(std::complex<double> ((mu(j) - lambda(i) + II * c)
          / (lambda(j) - lambda(i) + II * c)));
    }
    return result;
  }

}
