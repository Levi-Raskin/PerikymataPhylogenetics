#include "Eigen/Dense"
#include "Msg.hpp"
#include "ParameterInvVCV.hpp"
#include "ParameterCorrelationMatrix.hpp"
#include "ParameterDouble.hpp"
#include "ParameterRates.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TicToc.hpp"
#include "Utility.hpp"

#include <chrono>
#include <iostream>
#include <random>


ParameterInvVCV::ParameterInvVCV(double prob, std::string n, int numberOfTraits, PhylogeneticModel* p) : ParameterMatrix(prob, p, n), numtraits(numberOfTraits){
    
    //instantiate relative rates
    rates = new ParameterRates(1.0, model, "rates", numtraits);
    
    //instantiate correlation matrix
    correlationMatrix = new ParameterCorrelationMatrix(1.0, model, "corrmat", numtraits);
    components = {rates, correlationMatrix};
    double sum = 0.0;
    for (Parameter* p : components)
        sum += p->getProposalProbability();
    for (Parameter* p : components)
        p->setProposalProbability(p->getProposalProbability()/sum);
}

double ParameterInvVCV::getAcceptanceRatio(void){
    double sum = 0.0;
    for(Parameter* p : components)
        sum += p->getAcceptanceRatio();
    return sum / 3;
}

bool ParameterInvVCV::getAdaptiveProposalActive(void){
    for(Parameter* p : components)
        if(p->getAdaptiveProposalActive() == true)
            return true;
    return false;
}

Eigen::MatrixXd ParameterInvVCV::getValue(void){
    Eigen::VectorXd relRates = rates->getValue();
    Eigen::MatrixXd corr = correlationMatrix->getValue();
    Eigen::MatrixXd vcv = (relRates.asDiagonal() * corr * relRates.asDiagonal());
    
    return vcv.inverse();
}

double ParameterInvVCV::lnProbability(void){
    double lnPro = 0.0;
    
    for(Parameter* p : components)
        lnPro += p->lnProbability();
    
    Eigen::VectorXd relRates = rates->getValue();
    
    // Log Jacobian: n*log(avgRate) + 2*sum(log(relRates))
    double logJacobian = 2.0 * relRates.array().log().sum();
    
    return lnPro + logJacobian;
}

void ParameterInvVCV::print(void){
    std::cout << std::endl;
}

double ParameterInvVCV::update(void){
    Parameter* parm = nullptr;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double sum = 0.0;
    double u = rng.uniformRv();
    for (Parameter* p : components)
        {
        sum += p->getProposalProbability();
        if (u <= sum)
            {
            parm = p;
            break;
            }
        }
    updatedComponent = parm;
//    for(Parameter* p : components)
//        std::cout << p->getAcceptanceRatio() << " ";
//    std::cout << std::endl;
    return updatedComponent->update();
}

void ParameterInvVCV::updateForAcceptance(void){
    updatedComponent->updateForAcceptance();
}

void ParameterInvVCV::updateForRejection(void){
    updatedComponent->updateForRejection();
}

