#include "Msg.hpp"
#include "ParameterDouble.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

ParameterDouble::ParameterDouble(double prob, std::string n) : Parameter(prob, n), lowerBound(std::numeric_limits<double>::lowest() / 2), upperBound(std::numeric_limits<double>::max() / 2), numRejections(0), numAcceptances(0), numAdaptive(10000){
    adaptiveProposalActive = true;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double u = Probability::Normal::rv(&rng); //if no bounds provided, draw from standard normal
    value.push_back(u);
    value.push_back(u);
    windowSize = 1;
}

double ParameterDouble::lnProbability(void){
    double lnPro = -1.0;
    
    //Levi you should improve this
    if(lowerBound != 0){
        //placing a normal prior distribution if unbounded
        lnPro = Probability::Normal::lnPdf(0, 100.0, getValue());
    }else if(lowerBound == 0){
        //placing a gamma prior distirbution if bounded > 0
//        lnPro = Probability::Gamma::lnPdf(1.0, 1.0, getValue());
//        lnPro = 0.0;
    }else{
        //tbd
    }
    if(lnPro == -1.0)
        Msg::error("ParmDouble " + parmName + " prior probabiltiy is -1.0");
    
    return lnPro;
}

void ParameterDouble::print(void){

}

double ParameterDouble::update(void) {
    return updateAdaptive(10000, 0.43);
}

double ParameterDouble::updateAdaptive(int numGen, double targetR){
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    double acceptRej = 0.0;
    for(bool b : recentAcceptRej)
        if(b == true)
            acceptRej++;
    acceptRej /= recentAcceptRej.size();
    
    if((numRejections + numAcceptances) % 100 ==0 && ((numRejections + numAcceptances) < numGen)){
            if(acceptRej < targetR - 0.2)
                windowSize /= 1.1;
            else if (acceptRej > targetR + 0.2)
                windowSize *= 1.1;
            
    }else if ((numRejections + numAcceptances) == numGen){
        std::cout << parmName << " done adaptive sampling | final acceptRej: " << acceptRej << std::endl;
        adaptiveProposalActive = false;
    }

    double u = Probability::Uniform::rv(&rng, value[1] - windowSize, value[1] + windowSize);
    if(u <= lowerBound){ //bounce
        u = lowerBound + (lowerBound - u);
    }else if ( u >= upperBound){
        u = upperBound - ( u - upperBound);
    }
    value[0] = u;
    return 0.0;
}

void ParameterDouble::updateForAcceptance(void) {
    numAcceptances++;
    value[1] = value[0];
    recentAcceptRej.push_back(true);
    if(recentAcceptRej.size() > 1000)
            recentAcceptRej.pop_front();
}

void ParameterDouble::updateForRejection(void) {
    numRejections++;
    value[0] = value[1];
    recentAcceptRej.push_back(false);
    if(recentAcceptRej.size() > 1000)
            recentAcceptRej.pop_front();
}

void ParameterDouble::updateForRejectionNotDynamic(void) {
    value[0] = value[1];
}
