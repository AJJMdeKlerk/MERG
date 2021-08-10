#ifndef LiebLinState_H
#define LiebLinState_H

namespace ABACUS_VB {

  struct ScanInfo;

  class LiebLinState {
    private:
      SystemInfo systemInfo; 
      Eigen::VectorXi doubledQuantumNumbs;
      Eigen::VectorXd rapidities;
      Eigen::MatrixXd gaudinMatrix;
      double energy;
      double momentum;
      double thirdMoment;
      double lnNorm;

      double maxErrorRapidities;
      int maxNewtonIterations;
      bool allComputed;

      void checkPositiveInteractionStrength();
      void computeLnNorm();
      void computeMoments();
      void computeRapidities();
      void computeGaudinMatrix();
      void defaultInitRapidities();
      void defaultInitDoubledQuantumNumbs();
      void defaultInitAccuracy();
      Eigen::VectorXd BAEs(const Eigen::VectorXd& trialSolution);
      Eigen::VectorXd interactionDerivativeBAEs
        (const Eigen::VectorXd& trialSolution);
      void computeAll();
    public:
      void displayState() const;

      SystemInfo getSystemInfo() const;
      int getParticleNumber() const;
      double getInteractionStrength() const;
      double getSystemLength() const;
      double getMaxErrorRapidities() const;
      int getMaxNewtonIterations() const;
      int getSumDoubledQuantumNumbs() const;
      Eigen::VectorXi getDoubledQuantumNumbs() const;
      Eigen::VectorXd getRapidities() const;
      Eigen::MatrixXd getGaudinMatrix() const;
      SystemInfo getsystemInfo() const;
      double getRapidity(int i) const;
      double getEnergy() const;
      double getMomentum() const;
      double getThirdMoment() const;
      double getLnNorm() const;
      std::string getInfoForFilename() const;

      double kernel(const int i, const int j) const;

      void setDoubledQuantumNumbs(const Eigen::VectorXi& input);
      LiebLinState();
      LiebLinState(const SystemInfo& inputsystemInfo, const
          ScanInfo& ScanInfo, const std::string& label);
      LiebLinState(const SystemInfo& inputsystemInfo);
      LiebLinState(const SystemInfo& inputsystemInfo, const
          Eigen::VectorXi& inputDoubledQuantumNumbers);
  };

}

#endif
