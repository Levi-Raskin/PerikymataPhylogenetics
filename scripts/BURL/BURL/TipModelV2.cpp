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
    lnLDirty(true), lnPDirty(true),
    hasMissingData(false), gibbsPkUpdate(false){
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
    
    std::map<std::vector<int>, int> patternIdByKey;

    for(int i = 0; i < numRows; i++){
        if(tipDataIncomplete.row(i).array().isNaN().any()){
            std::vector<int> missIdx, obsIdx;
            missIdx.reserve(numCols);
            obsIdx.reserve(numCols);
            for(int j = 0; j < numCols; j++){
                if(std::isnan(tipDataIncomplete(i,j)))
                    missIdx.push_back(j);
                else
                    obsIdx.push_back(j);
            }

            auto it = patternIdByKey.find(missIdx);
            int patternIdx;
            if(it == patternIdByKey.end()){
                patternIdx = (int)patternMissingIdx.size();
                patternIdByKey.emplace(missIdx, patternIdx);
                patternMissingIdx.push_back(missIdx);
                patternObsIdx.push_back(obsIdx);
                patternRows.push_back({});
            }else{
                patternIdx = it->second;
            }
            patternRows[patternIdx].push_back(i);
        }
    }
    
    if(numRows == 1)
        Msg::warning("One observation given for " + tn + " | treating as species mean known without uncertainty");
    if(numRows == 1 && missingPkVals.empty() == false)
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
        
        updatableParameters.clear();
        for (Parameter* p : parameters)
            if (p->getProposalProbability() > 0.0)
                updatableParameters.push_back(p);
        
        //preallocs
        mu = taxonMean->getValue();
        sigma = taxonVariance->getValue();
        sigmaChol = sigma.llt();
        sigmaLogDet = sigmaChol.matrixLLT().diagonal().array().log().sum();\

        xDiff.resize(numCols);
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
    sigmaLogDet = sigmaChol.matrixLLT().diagonal().array().log().sum();
    
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
//        if(hasMissingData == true)
//            std::cout << " | missing data imputation a/r: " << (double)numImputationAcceptances / (double)(numImputationRejections + numImputationAcceptances);
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
    }else{
        gibbsPkUpdate = false;
        Parameter* parm = nullptr;
        double sum = 0.0;
        double u = rng.uniformRv();
        for (Parameter* p : updatableParameters)
            {
            sum += p->getProposalProbability();
            if (u < sum)
                {
                parm = p;
                break;
                }
            }
        if (parm == nullptr)
            parm = updatableParameters.back();
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

    if(patternRows.empty())
        Msg::error("updatePkGibbs called for " + tipName + " but no missingness patterns are present");

    // Pick a missingness-pattern group at random.
    int groupIdx = (int)(rng.uniformRv() * patternRows.size());
    if(groupIdx >= (int)patternRows.size())
        groupIdx = (int)patternRows.size() - 1;   // guard against uniformRv() == 1.0

    const std::vector<int>& missingIdx  = patternMissingIdx[groupIdx];
    const std::vector<int>& obsIdx      = patternObsIdx[groupIdx];
    const std::vector<int>& rowsInGroup = patternRows[groupIdx];

    const int numMissing = (int)missingIdx.size();
    const int numObs     = (int)obsIdx.size();

    if(obsIdx.empty())
        Msg::error("all tip data is missing for a row in " + tipName);

    sigma11.resize(numMissing, numMissing);
    sigma12.resize(numMissing, numObs);
    sigma21.resize(numObs, numMissing);
    sigma22.resize(numObs, numObs);
    u1.resize(numMissing);
    u2.resize(numObs);
    x2.resize(numObs);

    for(int i = 0; i < numMissing; i++){
        int mi = missingIdx[i];
        u1(i) = tipMean(mi);
        for(int j = 0; j < numMissing; j++)
            sigma11(i, j) = tipVCV(mi, missingIdx[j]);
        for(int j = 0; j < numObs; j++)
            sigma12(i, j) = tipVCV(mi, obsIdx[j]);
    }
    for(int i = 0; i < numObs; i++){
        int oi = obsIdx[i];
        u2(i) = tipMean(oi);
        for(int j = 0; j < numMissing; j++)
            sigma21(i, j) = tipVCV(oi, missingIdx[j]);
        for(int j = 0; j < numObs; j++)
            sigma22(i, j) = tipVCV(oi, obsIdx[j]);
    }

    // Factor sigma22 ONCE for this pattern group — reused for every row below.
    sigma22Solver = sigma22.llt();
    if(sigma22Solver.info() != Eigen::Success)
        Msg::error("Cholesky decomposition failed for " + tipName);

    sigma22Inv_sigma21 = sigma22Solver.solve(sigma21);
    sigmaCond = sigma11 - sigma12 * sigma22Inv_sigma21;

    for(int rowToUpdate : rowsInGroup){
        for(int i = 0; i < numObs; i++)
            x2(i) = tipDataIncomplete(rowToUpdate, obsIdx[i]);

        x2_minus_u2 = x2 - u2;
        sigma22Inv_diff = sigma22Solver.solve(x2_minus_u2);
        uCond = u1 + sigma12 * sigma22Inv_diff;

        newVals = Probability::MultivariateNormal::rv(&rng, uCond, &sigmaCond);

        for(int idx = 0; idx < numMissing; idx++){
            int j = missingIdx[idx];
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
