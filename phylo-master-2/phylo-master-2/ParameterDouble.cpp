#include "Msg.hpp"
#include "ParameterDouble.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include </usr/local/include/omp.h>

ParameterDouble::ParameterDouble(double prob, PhylogeneticModel* p,  std::string n, double min):  Parameter(prob, p, n), lowerBound(min), upperBound(std::numeric_limits<double>::max() / 2), numRejections(0), numAcceptances(0), numAdaptive(10000){
    adaptiveProposalActive = true;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double u = min + Probability::Exponential::rv(&rng, 1) + lowerBound;
    value.push_back(u);
    value.push_back(u);
    windowSize = 1;
}

ParameterDouble::ParameterDouble(double prob, PhylogeneticModel* p,  std::string n) : Parameter(prob, p, n), lowerBound(std::numeric_limits<double>::lowest() / 2), upperBound(std::numeric_limits<double>::max() / 2), numRejections(0), numAcceptances(0), numAdaptive(10000){
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
    if(lowerBound == 0)
        return updateExp(numAdaptive, 0.43);
//    return updateAdaptive(10000, 0.43);
    return updateNormal(numAdaptive, 0.43);
    /*
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    //sliding scale on value 0
    double acceptRej = ((double)numAcceptances)/((double)(numAcceptances+numRejections));
    if((numRejections + numAcceptances) % 5 ==0 && ((numRejections + numAcceptances) < 1000)){
            if((0.1 < acceptRej) && (acceptRej < 0.42))
                windowSize *= 0.8;
            else if (acceptRej < 0.1)
                windowSize *= 0.05;
            else if (acceptRej > 0.8)
                windowSize *= 2;
            else if (acceptRej > 0.44)
                windowSize *= 1.25;
            
    }else if ((numRejections + numAcceptances) == 1000){
        std::cout << parmName << " done adaptive sampling | final acceptRej: " << acceptRej << std::endl;
    }
    if(windowSize > 100)
        windowSize = 100; //window size is pegged at 100; seems to want a huge wide search area
//    if((numRejections + numAcceptances) % 5 ==0){
//        std::cout << std::setprecision(10);
//        std::cout  << "parm double sliding scale: " << windowSize << " | accept rej: " << acceptRej << std::endl;
//    }
     
    //windowsize is actually 1/2 true window size, this just saves division
    if(std::isnan(lowerBound) == true && std::isnan(upperBound) == true){
        value[0] = Probability::Uniform::rv(&rng, value[1] - windowSize, value[1] + windowSize);
    }else if(std::isnan(lowerBound) == false && std::isnan(upperBound) == true){
        double u = Probability::Uniform::rv(&rng, value[1] - windowSize, value[1] + windowSize);
        if(u <= lowerBound) //bounce
            u = lowerBound + (lowerBound - u);
        value[0] = u;
    }else{
        //tbd
    }
    //returns proposal ratio
    //for sliding scale, that is 1.0
    return 0.0;
    */
}

double ParameterDouble::updateAdaptive(int numGen, double targetR){
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    //sliding scale on value 0;

//    double acceptRej = ((double)numAcceptances)/((double)(numAcceptances+numRejections));
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

    //windowsize is actually 1/2 true window size, this just saves division
//    if(std::isnan(lowerBound) == true && std::isnan(upperBound) == true){
//        value[0] = Probability::Uniform::rv(&rng, value[1] - windowSize, value[1] + windowSize);
//    }else if(std::isnan(lowerBound) == false && std::isnan(upperBound) == true){
    double u = Probability::Uniform::rv(&rng, value[1] - windowSize, value[1] + windowSize);
    if(u <= lowerBound){ //bounce
        u = lowerBound + (lowerBound - u);
//        Msg::warning("we're up at the lower bound: "  + std::to_string(lowerBound));
    }else if ( u >= upperBound){
        u = upperBound - ( u - upperBound);
//        Msg::warning("we're up at the upper bound: " + std::to_string(upperBound));
    }
    value[0] = u;
//    }else{
//        //tbd
//        Msg::error("not implemented");
//    }
    //returns proposal ratio
    //for sliding scale, that is 1.0
    
//    if((numAcceptances + numRejections ) % 100 == 0)
//        std::cout <<  " recent accept rej: " << acceptRej << std::endl;
//    std::cout << windowSize << std::endl;
    return 0.0;
}

double ParameterDouble::updateExp(int numGen, double targetR){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    //sliding scale on value 0;

//    double acceptRej = ((double)numAcceptances)/((double)(numAcceptances+numRejections));
    double acceptRej = 0.0;
    for(bool b : recentAcceptRej)
        if(b == true)
            acceptRej++;
    acceptRej /= recentAcceptRej.size();
    
    if((numRejections + numAcceptances) % 10 ==0 && ((numRejections + numAcceptances) < numGen)){
            if(acceptRej < (targetR - 0.2))
                windowSize /= 1.1;
            else if (acceptRej > (targetR + 0.2))
                windowSize *= 1.1;
    }else if ((numRejections + numAcceptances) == numGen){
        std::cout << parmName << " done adaptive sampling | final acceptRej: " << acceptRej << std::endl;
        adaptiveProposalActive = false;
    }

    //windowsize is actually repurposed as lambda here
    value[0] = value[1] * std::exp(windowSize * (rng.uniformRv() - 0.5) );
    
    return std::log(value[0] / value[1]);
}

double ParameterDouble::updateNormal(int numGen, double targetR){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    //sliding scale on value 0;

//    double acceptRej = ((double)numAcceptances)/((double)(numAcceptances+numRejections));
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

    //windowsize is actually repurposed as lambda here
    value[0] = value[1] + Probability::Normal::rv(&rng, 0.0, windowSize);
    double n2o = Probability::Normal::lnPdf(value[0], windowSize, value[1]);
    double o2n = Probability::Normal::lnPdf(value[1], windowSize, value[0]);

    double hr =  n2o - o2n;
    
    return hr;
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
