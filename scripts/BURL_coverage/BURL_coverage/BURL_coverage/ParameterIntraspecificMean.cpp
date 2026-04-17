#include "ParameterIntraspecificMean.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include <iostream>
#include <set>

ParameterIntraspecificMean::ParameterIntraspecificMean(double prob, std::string n, Eigen::MatrixXd* data, PhylogeneticModel* p) : Parameter(prob, n), nTraits(data->cols()), nObs(data->rows()), tipData(data), windowSize(0.5), stepSize(1e-3), model(p), targetAcceptanceRate(0.43), upperAcceptanceRate(targetAcceptanceRate + 0.01),lowerAcceptanceRate(targetAcceptanceRate-0.01), ngAdaptive(50000), covarianceUpdateFreq(100), useEmpiricalCovariance(false), numAcceptances(0), numRejections(0){

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    //mean is unconstrained pos/neg
    //prior
    mu0 = Eigen::VectorXd::Ones(nTraits);
    for(int i = 0; i < nTraits; i++)
        mu0(i) = Probability::Normal::rv(&rng, 1.0, 1.0);

    mean.push_back(tipData->colwise().mean());
    mean.push_back(tipData->colwise().mean());
    adaptiveProposalActive = true;
    
    dof = nTraits + 2; //such that the mean is the scale matrix
    lambda = dof - nTraits - 1;

    mu0 = Eigen::VectorXd::Zero(nTraits);
    
    psi = Eigen::MatrixXd::Constant(nTraits, nTraits, 1e-6);
    psi.diagonal().array() = 1.0;
    psi = Eigen::VectorXd::Ones(nTraits).asDiagonal() * psi * Eigen::VectorXd::Ones(nTraits).asDiagonal();
    
    
    empiricalMean = Eigen::VectorXd::Zero(nTraits);
    empiricalCovariance = Eigen::MatrixXd::Identity(nTraits, nTraits);
    
    proposalCov.resize(nTraits, nTraits);
}

double ParameterIntraspecificMean::lnProbability(void){
    //not needed; accounted for by multivariateBrownain Motion
    return 0.0;
}

void ParameterIntraspecificMean::print(void){
    std::cout << " " << std::endl;
}

double ParameterIntraspecificMean::update(void){
    return updateMHmvNDraw();
}

double ParameterIntraspecificMean::updateMHmvNDraw(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    const size_t totalSamples = numRejections + numAcceptances;
    
    if (totalSamples % 50 == 0 && totalSamples < ngAdaptive) {
    const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
        const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
        if (acceptanceRate < lowerAcceptanceRate)
            windowSize /= 1.05;
        else if (acceptanceRate > upperAcceptanceRate)
            windowSize *= 1.05;
//        std::cout << windowSize << "\n";
    }else if ((numRejections + numAcceptances) == ngAdaptive){
        adaptiveProposalActive = false;
    }
    
     if (totalSamples % covarianceUpdateFreq == 0 && totalSamples < ngAdaptive) {
        updateEmpiricalCovariance();
    }
    
    if (totalSamples == 5000 && !useEmpiricalCovariance) {
        useEmpiricalCovariance = true;
    }
    
    // Choose proposal covariance
    if(totalSamples < ngAdaptive){
        if (useEmpiricalCovariance && recentSamples.size() > nTraits * 2) {
            proposalCov = windowSize * empiricalCovariance;
            proposalCholLower = proposalCov.llt().matrixL().toDenseMatrix();
        } else {
            proposalCov = windowSize * Eigen::MatrixXd::Identity(nTraits, nTraits);
            proposalCholLower = proposalCov.llt().matrixL().toDenseMatrix();
        }
    }
    
    // Propose
    Probability::MultivariateNormal::rv(&rng, mean[0], mean[1], proposalCholLower);
    return 0.0;

}

void ParameterIntraspecificMean::updateEmpiricalCovariance(void){
    if (recentSamples.size() < nTraits + 1) return;
    
    // Compute empirical mean
    empiricalMean = Eigen::VectorXd::Zero(nTraits);
    for (const auto& sample : recentSamples) {
        empiricalMean += sample;
    }
    empiricalMean /= recentSamples.size();
    
    // Compute empirical covariance
    empiricalCovariance = Eigen::MatrixXd::Zero(nTraits, nTraits);
    for (const auto& sample : recentSamples) {
        Eigen::VectorXd centered = sample - empiricalMean;
        empiricalCovariance += centered * centered.transpose();
    }
    empiricalCovariance /= (recentSamples.size() - 1);
    
    // Add small diagonal for numerical stability
    empiricalCovariance += 1e-6 * Eigen::MatrixXd::Identity(nTraits, nTraits);
}

void ParameterIntraspecificMean::updateForAcceptance(void){
    mean[1] = mean[0];
    numAcceptances++;
    recentAcceptRej.push_back(true);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
    recentSamples.push_back(mean[0]);
    if (recentSamples.size() > 1000)  // Keep last 500 samples
        recentSamples.pop_front();
}

void ParameterIntraspecificMean::updateForRejection(void){
    mean[0] = mean[1];
    numRejections++;
    recentAcceptRej.push_back(false);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}
