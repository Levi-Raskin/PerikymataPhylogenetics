#include "Msg.hpp"
#include "ParameterRates.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Utility.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include </usr/local/include/omp.h>

ParameterRates::ParameterRates(double prob, PhylogeneticModel* p, std::string n, int nt) : Parameter(prob, p, n), numRejections(0), numAcceptances(0), windowSize(0.1), nTraits(nt){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Eigen::VectorXd logRates = Eigen::VectorXd::Zero(nTraits);
    for(int i = 0; i < nTraits; i++)
        logRates(i) = Probability::Normal::rv(&rng);
    value.push_back(logRates) ;
    value.push_back(logRates) ;
}

double ParameterRates::lnProbability(void){
    double pr = 0.0;
    for(int i = 0; i < nTraits; i++)
        pr += Probability::Normal::lnPdf(0.0, 1.0, value[0](i));
    return pr;
}

void ParameterRates::print(void){
    Utility::EigenUtils::printEigen(value[0].array().exp());
}

double ParameterRates::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    if(rng.uniformRv() < 0.25){
        return updateHMC();
    }else{
        return updateMH();
    }
}

double ParameterRates::updateMH(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    double acceptRej = 0.0;
    for(bool b : recentAcceptRej)
        if(b == true)
            acceptRej++;
    acceptRej /= recentAcceptRej.size();
    
    if((numRejections + numAcceptances) % 100 ==0 && ((numRejections + numAcceptances) < 10000)){
            if(acceptRej < 0.48)
                windowSize /= 1.1;
            else if (acceptRej > 0.38)
                windowSize *= 1.1;
    }

    //sample a random element to update
    int element = (int)(rng.uniformRv() * nTraits);

    double u = Probability::Uniform::rv(&rng, value[0](element) - windowSize, value[0](element) + windowSize);
    value[0](element) = u;
    
    return 0.0;
}

double ParameterRates::updateHMC(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
//    t.end();

    // ---- 1. Flatten parameters: rates + Cholesky off-diagonals ----
    Eigen::VectorXd flatParams(nTraits);
    flatParams << value[0];

    // ---- 2. Define log probability ----
    auto logProb = [&](const Eigen::VectorXd& theta) -> double {
        double lp = model->lnPriorProbability();
        lp += model->lnLikelihood();
        return lp;
    };

    // ---- 3. Numerical gradient ----
    auto numericalGradient = [&](const Eigen::VectorXd& theta, double eps = 1e-6) {
        Eigen::VectorXd grad(nTraits);
//        #pragma omp parallel for num_threads(10)
            for (int i = 0; i < nTraits; ++i) {
                Eigen::VectorXd perturbed = theta;
                perturbed(i) += eps;
                double f1 = logProb(perturbed);
                perturbed(i) -= 2 * eps;
                double f2 = logProb(perturbed);
                grad(i) = (f1 - f2) / (2 * eps);
            }
        return grad;
    };

    // ---- 4. Draw momenta ----
    Eigen::VectorXd momentum = Eigen::VectorXd::NullaryExpr(nTraits, [&]() {
        return Probability::Normal::rv(&rng);
    });

    // ---- 5. Compute initial energy ----
    Eigen::VectorXd gradU0 = -numericalGradient(flatParams);

    // ---- 6. Leapfrog integration ----
    double eps = 1e-3;
        
    int L = 10;
    Eigen::VectorXd theta = flatParams;
    Eigen::VectorXd p = momentum;

    p -= 0.5 * eps * gradU0;
    for (int i = 0; i < L; ++i) {
        theta += eps * p;
        value[0] = theta;
        Eigen::VectorXd gradU = -numericalGradient(theta);
        if (i != L - 1)
            p -= eps * gradU;
        else
            p -= 0.5 * eps * gradU;
    }
//    t.end();
    return 0.0;
}

void ParameterRates::updateForAcceptance(void){
    numAcceptances++;
    value[1] = value[0];
    recentAcceptRej.push_back(true);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}

void ParameterRates::updateForRejection(void){
    numRejections++;
    value[0] = value[1];
    recentAcceptRej.push_back(false);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}
