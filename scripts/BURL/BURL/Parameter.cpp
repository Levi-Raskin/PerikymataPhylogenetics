#include "Parameter.hpp"



Parameter::Parameter(double prob, PhylogeneticModel* m, std::string n) : proposalProbability(prob), model(m), parmName(n), adaptiveProposalActive(false), parmPrintsToConsole(true) {

}

Parameter::Parameter(double prob, std::string n) : proposalProbability(prob), parmName(n), adaptiveProposalActive(false), parmPrintsToConsole(true) {

}
