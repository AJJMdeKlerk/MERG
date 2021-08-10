#ifndef MatrixElements_H
#define MatrixElements_H

namespace ABACUS_VB {

  double rho(const LiebLinState& leftState, const LiebLinState& rightState);
  std::complex<double> lnRhoOffDiagonal(const LiebLinState& leftState, 
    const LiebLinState& rightState); 
  double G2(const LiebLinState& leftState, const LiebLinState& rightState);
  std::complex<double> lnG2OffDiagonal(const LiebLinState& leftState, 
    const LiebLinState& rightState); 
  std::complex<double> lnVPlus(const Eigen::VectorXd& mu, 
      const Eigen::VectorXd& lambda, double c, int i);
  std::complex<double> lnVPlusArbitrary(const Eigen::VectorXd& mu, 
      const Eigen::VectorXd& lambda, std::complex<double> lambdaP,
      double c);

  double G2Diagonal(const LiebLinState& state);
  Eigen::VectorXd computeInteractionDerivativeRapidities 
    (const LiebLinState& state);

  Eigen::VectorXd initInteractionDerivativeRapidities 
    (const LiebLinState& state);
  Eigen::VectorXd interactionDerivativeBAEs 
    (const LiebLinState& state, const Eigen::VectorXd& trialSolution);

  bool isCompatible(const LiebLinState& leftState, 
    const LiebLinState& rightState);
 
  double kernel(const double& lambda, const double& mu, const double& c);

  template <typename T>
  std::complex<double> lnDetInvertibleMatrix(T input) {
    Eigen::PartialPivLU<T> luDecomp(input);
    Eigen::MatrixXd permutation = luDecomp.permutationP();
    std::complex<double> lnNorm = 0.0;
    lnNorm += std::log(std::complex<double>(permutation.determinant()));
    for (int i = 0; i < input.rows(); i++) {
      lnNorm += std::log(std::complex<double>(luDecomp.matrixLU()(i, i)));
    }
    return lnNorm;
  }

}

#endif




