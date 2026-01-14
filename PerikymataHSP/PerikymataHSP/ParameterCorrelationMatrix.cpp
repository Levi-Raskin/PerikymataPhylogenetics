#include "Msg.hpp"
#include "ParameterCorrelationMatrix.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Utility.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include </usr/local/include/omp.h>

ParameterCorrelationMatrix::ParameterCorrelationMatrix(double prob, PhylogeneticModel* p, std::string n, int nt) : Parameter(prob, p, n), numRejections(0), numAcceptances(0), nTraits(nt), windowSize(1){

    Eigen::MatrixXd correlationMatrix0 = Eigen::MatrixXd::Constant(nTraits, nTraits, 1e-6);
    correlationMatrix0.diagonal().array() = 1.0;
//    Eigen::LLT<Eigen::MatrixXd> llt = correlationMatrix0.llt();
//    Eigen::MatrixXd L = llt.matrixL();
//
//    value.push_back(L);
//    value.push_back(L);
    value.push_back(correlationMatrix0);
    value.push_back(correlationMatrix0);
}

double ParameterCorrelationMatrix::lnProbability(void){
    return 0.0;
}

void ParameterCorrelationMatrix::print(void){
    Utility::EigenUtils::printEigen(value[0]);
}

double ParameterCorrelationMatrix::update(void){
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
    int col = (int)(rng.uniformRv() * nTraits);
    int row = (int)(rng.uniformRv() * nTraits);
    
    while(row == col)
        col = (int)(rng.uniformRv() * nTraits);
    
    double current_value = value[0](row,col);
    double u  = Probability::TruncatedNormal::rv(&rng, current_value, windowSize, -1.0, 1.0);

    value[0](row,col) = u;
    value[0](col,row) = u;
    
    double hr = Probability::TruncatedNormal::lnPdf(current_value, u, windowSize, -1.0, 1.0) - Probability::TruncatedNormal::lnPdf(u, current_value, windowSize, -1.0, 1.0); // backwards - forwards
    return hr;
    
    // transform the current value from [-1, 1] to [0, 1]
//    current_value = (current_value + 1.0) / 2.0;
//    
//    double alpha = windowSize;
//    
//    // draw new rates and compute the hastings ratio at the same time
//    double a = alpha + 1.0;
//    double b = (a - 1.0) / current_value - a + 2.0;
//
//    double new_value = Probability::Beta::rv(&rng, a, b);
//
//    // set the value (for both sides of the matrix!)
//    double new_value_transformed = new_value * 2.0 - 1.0;
//       
//    value[0](row,col) = new_value_transformed;
//    value[0](col,row) = new_value_transformed;
//    
//    double ln_Hastings_ratio = 0.0;
//    
//    // compute the Hastings ratio
//    double forward = Probability::Beta::lnPdf(a, b, new_value);
//    double new_a = alpha + 1.0;
//    double new_b = (a - 1.0) / new_value - a + 2.0;
//    double backward = Probability::Beta::lnPdf(new_a, new_b, current_value);
//    ln_Hastings_ratio = backward - forward;
//    return ln_Hastings_ratio;
}

void ParameterCorrelationMatrix::updateForAcceptance(void){
    numAcceptances++;
    value[1] = value[0];
    recentAcceptRej.push_back(true);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}

void ParameterCorrelationMatrix::updateForRejection(void){
    numRejections++;
    value[0] = value[1];
    recentAcceptRej.push_back(false);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}
