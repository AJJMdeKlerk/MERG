#include "ABACUS_VB.h"

using namespace ABACUS_VB;

int main() {

  // Set the parameters of the Hamiltonian after the quench
  int particleNumber = 10;
  double systemLength = 10.0;
  double interactionStrength = 10.0;
  SystemInfo systemInfo {particleNumber, interactionStrength, systemLength};

  // Initialize a Bethe state, in this case the ground state
  LiebLinState groundState(systemInfo);

  // To consider an excited state, initialize a state as 
  // LiebLinState excitedState (systemInfo, doubledQuantumNumbs);
  // where doubledQuantumNumbs is an object of the type Eigen::VectorXi
  // specifying the quantum numbers multiplied by two

  // Set the seed state for the scanning routine, similarly the state that the
  // MERG-routine will track if that option is chosen later on
  LiebLinState seedState = groundState;
  std::cout << "The chosen seed state is: " << std::endl;
  seedState.displayState();

  // Set the parameters for the NRG-routine
  int statesKept = 75;
  int statesAdded = 25;
  int maxNRGSteps = 100;
  std::string method = "MERG";
  NRGInfo nrgInfo = {statesKept, statesAdded, maxNRGSteps, method};

  // Set the parameters for the state we are constructing in the basis of the
  // Hamiltonian after the quench. In this case, the ground state of the
  // Lieb-Liniger model at a different interaction strength
  SystemInfo systemInfoTargetState = {particleNumber, interactionStrength +
    changeInteractionStrength, systemLength};
  LiebLinState targetState(systemInfoTargetState);
  InfoTargetState infoTargetState = {targetState.getEnergy()};
  std::cout << "The target energy is " << targetState.getEnergy() << std::endl;
  std::cout << std::endl;
  std::string quench = "InteractionStrength";
  QuenchInfo quenchInfo = {quench, perturbationStrength, infoTargetState};

  // Set the parameters for the time evolution of the operator we want to
  // consider
  double totalTime = 1.0;
  double timeStepSize = 0.0001;
  std::string timeEvolvedOperator = "G2";
  TimeEvolutionInfo timeEvolutionInfo(totalTime, timeStepSize, timeEvolvedOperator);

  // Start the NRG-routine, which includes a call to the relevant preferential
  // scanning routine, see src/ConstructorNRG.cpp for an overview of the
  // algorithm
  NRG interactionQuench(quenchInfo, scanInfo, nrgInfo, timeEvolutionInfo);k

  return 0;
}
