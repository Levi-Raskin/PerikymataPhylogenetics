#include "Msg.hpp"
#include "ParameterDouble.hpp"
#include "RandomVariable.hpp"
#include "ParameterIntraspecificMean.hpp"
#include "ParameterIntraspecificVarianceGibbs.hpp"
#include "PerikymataHSPv4.hpp"
#include "Probability.hpp"
#include "TipModelV2.hpp"

#include <iostream>

TipModelV2::TipModelV2(std::string tn, Eigen::MatrixXd d, PerikymataHSPv4* m) : PhylogeneticModel(),
    tipDataIncomplete(d),
    tipName(tn),
    model(m),
    numCols((int)tipDataIncomplete.cols()),
    numRows((int)tipDataIncomplete.rows()),
    numImputationRejections(0), numImputationAcceptances(0),
    cachedLnL(0.0), cachedLnP(0.0),
    lnLDirty(true), lnPDirty(true){
    //Data wrangling and setting up tipdatacomplete and missingPkVals objects
    for(int i = 0; i < tipDataIncomplete.rows(); i++){
        for(int j = 0; j < tipDataIncomplete.cols(); j++){
            if(std::isnan(tipDataIncomplete(i,j))){
                ParameterDouble* newParm = new ParameterDouble(
                    0.0, 
                    "missing_" + tipName + "_(" + std::to_string(i) + "," + std::to_string(j) + ")"
                );
                newParm->setParmPrintConsole(false);
                parameters.push_back(newParm);
                missingPkVals.emplace(std::make_pair(i, j), newParm);
            }
        }
    }
    updatedImpPkDoubles.reserve(parameters.size()); // avoiding memory shuffling overhead later
    
    if(numRows == 1)
        Msg::warning("One observation given for " + tn + " | treating as species mean known without uncertainty");
    if(numRows == 1 && updatedImpPkDoubles.size() != 0)
        Msg::error("Imputation of species mean missing data not yet supported; coming soon");
    
    updateTipDataComplete();
    
    if(numRows > 1){
        taxonVariance = new ParameterIntraspecificVarianceGibbs(1.0, tipName + "_vcv", &tipDataComplete, this);
        taxonVariance->setParmPrintConsole(false);
        taxonMean = new ParameterIntraspecificMean(10.0, tipName + "_mean", &tipDataComplete, m);
        parameters.push_back(taxonMean);
        parameters.push_back(taxonVariance);
        taxonMean->setVarianceCovarianceMatrix(taxonVariance);
        taxonVariance->setMean(taxonMean);
        
        taxonVariance->update();
        taxonVariance->updateForAcceptance();
        
        hasMissingData = tipDataIncomplete.array().isNaN().any();
        
        double sum = 0.0;
        for (Parameter* p : parameters)
            sum += p->getProposalProbability();
        for (Parameter* p : parameters)
            p->setProposalProbability(p->getProposalProbability()/sum);
        
        
        //preallocs
        mu = taxonMean->getValue();
        sigma = taxonVariance->getValue();
        sigmaChol = sigma.llt();
        L = sigmaChol.matrixL();
        sigmaLogDet = L.diagonal().array().log().sum();
        
        xDiff.resize(numRows);
        scratchVec.resize(taxonMean->getValue().size());
        scratchMat.resize(taxonVariance->getValue().rows(), taxonVariance->getValue().cols());
        
        log2pi = std::log(2 * PI);
        numCnumR = numCols * numRows;
        term1 = -numCnumR / 2 * log2pi;
    }else{
        mu = tipDataIncomplete.row(0);
    }
}

TipModelV2::~TipModelV2(void){
    if(numRows > 1){
        for(auto& s : missingPkVals)
            delete s.second;
        delete taxonMean;
        delete taxonVariance;
    }
}

const Eigen::VectorXd& TipModelV2::getTipMean(void){
    if(numRows == 1)
        return mu;
    return taxonMean->getValue();
}

std::vector<std::string> TipModelV2::getParameterNames(void){
    parmNames.clear();
    for(int p = 0; p < parameters.size(); p++){
        ParameterIntraspecificMean* pim = dynamic_cast<ParameterIntraspecificMean*>(parameters[p]);
        ParameterIntraspecificVarianceGibbs* piv = dynamic_cast<ParameterIntraspecificVarianceGibbs*>(parameters[p]);
        if(pim != nullptr){
            for(int i = 0; i < numCols; i++)
                parmNames.push_back(pim->getName()+ "_" + std::to_string(i));
        }else if(piv != nullptr){
            for(int i = 0; i < numCols; i++)
                for(int j = 0; j < numCols; j++)
                    parmNames.push_back(piv->getName() + "_(" + std::to_string(i) + "," +  std::to_string(j) + ")");
        }else{
            parmNames.push_back(parameters[p]->getName());
        }
    }
    return parmNames;
}

std::vector<double> TipModelV2::getParameterString(void){
    parmValues.clear();
    for(int i = 0; i < parameters.size(); i++){
        ParameterDouble* pt = dynamic_cast<ParameterDouble*>(parameters[i]);
        ParameterIntraspecificMean* pim = dynamic_cast<ParameterIntraspecificMean*>(parameters[i]);
        ParameterIntraspecificVarianceGibbs* piv = dynamic_cast<ParameterIntraspecificVarianceGibbs*>(parameters[i]);
        if(pt != nullptr){
            parmValues.push_back(pt->getValue());
        }else if(pim != nullptr){
            scratchVec = pim->getValue();
            for(int i = 0; i < scratchVec.size(); i++)
                parmValues.push_back(scratchVec(i));
        }else if(piv != nullptr){
            scratchMat = piv->getValue();
            for(int i = 0; i < scratchMat.rows(); i++)
                for(int j = 0; j < scratchMat.cols(); j++)
                    parmValues.push_back(scratchMat(i,j));
        }else{
            parmValues.push_back(-1.0);
        }
    }
    return parmValues;
}

double TipModelV2::computeLnLikelihood(void){
    if(numRows == 1)
        return 0.0;
    mu = taxonMean->getValue();
    sigma = taxonVariance->getValue();
    sigmaChol = sigma.llt();
    L = sigmaChol.matrixL();
    sigmaLogDet = L.diagonal().array().log().sum();
    
    double term2 = -numRows * sigmaLogDet;
    
    double term3 = 0.0;
    for(int i = 0; i < numRows; i++){
        xDiff = tipDataComplete.row(i) - mu.transpose();
        term3 += xDiff.transpose() * sigmaChol.solve(xDiff);
    }
    term3 *= -0.5;
    double lnl = term1 + term2 + term3;
    return lnl;
}

double TipModelV2::computeLnPriorProbability(void){
    if(numRows == 1)
        return 0.0;
//    double lnp = 0.0;
//    for(auto p : parameters)
//        lnp += p->lnProbability();
    return taxonVariance->lnProbability();
}

double TipModelV2::lnLikelihood(void){
    if(lnLDirty){
        cachedLnL = computeLnLikelihood();
        lnLDirty  = false;
    }
    return cachedLnL;
}

double TipModelV2::lnPriorProbability(void){
    if(lnPDirty){
        cachedLnP = computeLnPriorProbability();
        lnPDirty  = false;
    }
    return cachedLnP;
}

void TipModelV2::print(void){
    if(numRows > 1){
        std::cout << " -- ";
        for(Parameter* p : parameters)
            if(p->getParmPrintConsole() == true)
                std::cout << p->getName() << " a/r: " << p->getAcceptanceRatio() << " " << p->getAdaptiveProposalActive();
        if(hasMissingData == true)
            std::cout << " | missing data imputation a/r: " << (double)numImputationAcceptances / (double)(numImputationRejections + numImputationAcceptances);
        std::cout << "\n";
    }
}

double TipModelV2::update(void){
    if(numRows == 1)
        return 0.0;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    if(hasMissingData == true && rng.uniformRv() < 0.1){
        gibbsPkUpdate = true;
        updatePkGibbs();
        lnLDirty = true;
        lnPDirty = true;
    return std::numeric_limits<double>::max();
        return 0.0;
    }else{
        gibbsPkUpdate = false;
        Parameter* parm = nullptr;
        double sum = 0.0;
        double u = rng.uniformRv();
        for (Parameter* p : parameters)
            {
            sum += p->getProposalProbability();
            if (u <= sum)
                {
                parm = p;
                break;
                }
            }
        updatedParameter = parm;
        double hr = updatedParameter->update();
        lnLDirty = true;
        lnPDirty = true;
        return hr;
    }
}

void TipModelV2::updateForAcceptance(void){
    if(numRows > 1){
        if(gibbsPkUpdate == true){
            //tipDataComplete is already updated; no need to calculate tip data again
            numImputationAcceptances++;
            for(ParameterDouble* p : updatedImpPkDoubles)
                p->updateForAcceptance();
            updateTipDataComplete();
        }else{
            updatedParameter->updateForAcceptance();
        }
    }
}

void TipModelV2::updateForRejection(void){
    if(numRows > 1){
        if(gibbsPkUpdate == true){
            numImputationRejections++;
            for(ParameterDouble* p : updatedImpPkDoubles)
                p->updateForRejection();
            updateTipDataComplete();
            lnLDirty = true;
            lnPDirty = true;
        }else{
            updatedParameter->updateForRejection();
            lnLDirty = true;
            lnPDirty = true;
        }
    }
}

void TipModelV2::updateTipDataComplete(void){
    tipDataComplete = tipDataIncomplete; //populated tip data (replaces NAN with associated parameter double value
    for(int i = 0; i < tipDataComplete.rows(); i++)
        for(int j = 0; j < tipDataComplete.cols(); j++)
            if(std::isnan(tipDataIncomplete(i,j))){
                auto key = std::make_pair(i, j);
                auto it = missingPkVals.find(key);
                if (it == missingPkVals.end() || it->second == nullptr)
                    Msg::error("missingPkVals missing key");
                tipDataComplete(i, j) = it->second->getValue();
            }
}

void TipModelV2::updatePkGibbs(void){
    const Eigen::VectorXd& tipMean = taxonMean->getValue();
    const Eigen::MatrixXd& tipVCV = taxonVariance->getValue();
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    updatedImpPkDoubles.clear();

    //sample a row to update (just updating one row at a time; less aggressive)
    int rowToUpdate = (int)(rng.uniformRv() * tipDataIncomplete.rows());
    
    const Eigen::VectorXd& indDat = tipDataIncomplete.row(rowToUpdate);

    // Count missing values first
    missingIndices.clear();
    obsIndices.clear();
    missingIndices.reserve(indDat.size());
    obsIndices.reserve(indDat.size());
    
    for(int i = 0; i < indDat.size(); i++) {
        if(std::isnan(indDat(i)))
            missingIndices.push_back(i);
        else
            obsIndices.push_back(i);
    }
    
    // Skip individuals with no missing data
    if(!missingIndices.empty()){
        if(obsIndices.empty())
            Msg::error("all tip data is missing");
        
        int numMissing = (int)missingIndices.size();
        int numObs = (int)obsIndices.size();
        
        // Resize matrices only when necessary (reuse memory)
        sigma11.resize(numMissing, numMissing);
        sigma12.resize(numMissing, numObs);
        sigma21.resize(numObs, numMissing);
        sigma22.resize(numObs, numObs);
        u1.resize(numMissing);
        u2.resize(numObs);
        x2.resize(numObs);
        
        for(int i = 0; i < numMissing; i++) {
            int mi = missingIndices[i];
            u1(i) = tipMean(mi);
            for(int j = 0; j < numMissing; j++) {
                sigma11(i, j) = tipVCV(mi, missingIndices[j]);
            }
            for(int j = 0; j < numObs; j++) {
                sigma12(i, j) = tipVCV(mi, obsIndices[j]);
            }
        }
        
        for(int i = 0; i < numObs; i++) {
            int oi = obsIndices[i];
            u2(i) = tipMean(oi);
            x2(i) = indDat(oi);
            for(int j = 0; j < numMissing; j++) {
                sigma21(i, j) = tipVCV(oi, missingIndices[j]);
            }
            for(int j = 0; j < numObs; j++) {
                sigma22(i, j) = tipVCV(oi, obsIndices[j]);
            }
        }
        
        sigma22Solver = sigma22.llt();
        if(sigma22Solver.info() != Eigen::Success) {
            Msg::error("Cholesky decomposition failed");
        }
        
        // Compute conditional distribution parameters
        sigma22Inv_sigma21 = sigma22Solver.solve(sigma21);
        sigmaCond = sigma11 - sigma12 * sigma22Inv_sigma21;
        
        // sigma22Inv * (x2 - u2)
        x2_minus_u2 = x2 - u2;
        sigma22Inv_diff = sigma22Solver.solve(x2_minus_u2);
        uCond = u1 + sigma12 * sigma22Inv_diff;
        
        // Sample new values for missing data
        newVals = Probability::MultivariateNormal::rv(&rng, uCond, &sigmaCond);
        for(int idx = 0; idx < numMissing; idx++) {
            int j = missingIndices[idx];
            auto key = std::make_pair(rowToUpdate, j);
            auto it = missingPkVals.find(key);
            if (it != missingPkVals.end() && it->second != nullptr) {
                it->second->setValue(newVals(idx));
                updatedImpPkDoubles.push_back(it->second);
            } else {
                Msg::error("Error: missing imputedPkVals entry");
            }
            tipDataComplete(rowToUpdate, j) = newVals(idx);
        }
    }
}
