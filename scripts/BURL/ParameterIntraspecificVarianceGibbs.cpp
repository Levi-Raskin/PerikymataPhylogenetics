#include "ParameterIntraspecificVarianceGibbs.hpp"
#include "ParameterIntraspecificMean.hpp"
#include "Eigen/Dense"
#include "Msg.hpp"
#include "ParameterVectorDouble.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TicToc.hpp"
#include "Utility.hpp"

#include <iostream>

ParameterIntraspecificVarianceGibbs::ParameterIntraspecificVarianceGibbs(double prob, std::string n, Eigen::MatrixXd* data, PhylogeneticModel* p) : ParameterMatrix(prob, n), numRejections(0), numAcceptances(0), numtraits(data->cols()), nObs(data->rows()), model(p), tipData(data), useCachedLnP(false), cachedlnP(-std::numeric_limits<double>::infinity()){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    dof = numtraits + 2; //such that the mean is the scale matrix
    lambda = dof - numtraits - 1;
    mu0 = Eigen::VectorXd::Zero(numtraits);
    psi = Eigen::MatrixXd::Constant(numtraits, numtraits, 1e-6);
    psi.diagonal().array() = 1.0;
    psi = Eigen::VectorXd::Ones(numtraits).asDiagonal() * psi * Eigen::VectorXd::Ones(numtraits).asDiagonal();
    
    value.push_back(Probability::InverseWishart::rv(&rng, psi, dof));
    value.push_back(value[0]);

}

bool ParameterIntraspecificVarianceGibbs::getAdaptiveProposalActive(void){
    return false;
}

const Eigen::MatrixXd& ParameterIntraspecificVarianceGibbs::getValue(void) {
    return value[0];
}

double ParameterIntraspecificVarianceGibbs::lnProbability(void){
    //prior probability of sigma
    if(useCachedLnP == false){
        cachedlnP = Probability::InverseWishart::lnPdf(&(value[0]), &psi, dof);
        useCachedLnP = true;
    }
    return cachedlnP;
}

void ParameterIntraspecificVarianceGibbs::print(void){
}

double ParameterIntraspecificVarianceGibbs::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Eigen::VectorXd mu_current = mean->getValue();
    
    double dofN = dof + nObs;
    
    //computing scatter aroudn current mean value conditions on m
    Eigen::MatrixXd datScatter = Eigen::MatrixXd::Zero(numtraits, numtraits);
    for(int i = 0; i < nObs; i++){
        Eigen::VectorXd ydiff = tipData->row(i).transpose() - mu_current;
        datScatter += ydiff * ydiff.transpose();
    }
    
    Eigen::VectorXd diff = mu_current - mu0;
    Eigen::MatrixXd prior_scatter = lambda * (diff * diff.transpose());
    Eigen::MatrixXd psiN = psi + datScatter + prior_scatter;
    Eigen::MatrixXd psiNInvLower = psiN.inverse().llt().matrixL();
    
//    value[0] = Probability::InverseWishart::rv(&rng, &psiN, dofN);
    Probability::InverseWishart::rv(&rng, value[0], psiNInvLower, dofN);
    useCachedLnP = false;
    return std::numeric_limits<double>::max(); //this is a gibbs update, hence should always be accepted
}

void ParameterIntraspecificVarianceGibbs::updateForAcceptance(void){
    numAcceptances++;
    value[1] = value[0];
}

void ParameterIntraspecificVarianceGibbs::updateForRejection(void){
    numRejections++;
    value[0] = value[1];
    useCachedLnP = false;
}
