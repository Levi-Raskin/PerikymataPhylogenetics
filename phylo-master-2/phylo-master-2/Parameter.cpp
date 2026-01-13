#include "Parameter.hpp"

Parameter::Parameter(double proposalProb, PhylogeneticModel* m, std::string n) : proposalProbability(proposalProb), model(m), parmName(n), adaptiveProposalActive(false) {

}
