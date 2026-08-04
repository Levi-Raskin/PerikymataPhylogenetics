#include "Msg.hpp"
#include "MultivariateBrownianMotionV2.hpp"
#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterIntraspecificMean.hpp"
#include "ParameterTree.hpp"
#include "PerikymataHSPv4.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TipModelV2.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <iostream>
#include <string>
#include <vector>

PerikymataHSPv4::PerikymataHSPv4(Tree* backbone, std::vector<std::string> datRN, Eigen::MatrixXd* dat) : MultivariateBrownianMotionV2(),
    fixedTree(*backbone),
    updateTipsOn(true),
    savedTipIdx(-1){
    
    UserSettings& settings = UserSettings::userSettings();
    if(settings.getWithIntraspecific() == false)
        updateTipsOn = false;
    withPhylo = settings.getWithPhylogeny();
    
    fixedTree = Tree(*backbone);

    // First pass: identify tips to drop
    std::vector<std::string> toDrop;
    for(Node* n : fixedTree.getDownPassSequence()){
        if(n->getIsTip()){
            bool found = false;
            for(const std::string& s : datRN){
                if(s == n->getName()){
                    found = true;
                    break;
                }
            }
            if(!found){
                toDrop.push_back(n->getName());
            }
        }
    }

    for(const std::string& name : toDrop){
        std::cout << "dropping " << name << std::endl;
        fixedTree.dropTip(name);
    }

    fixedTree.reindexNodes();
    fixedTree.initializeDownPassSequence();

    std::unordered_map<std::string, int> nodeIdxByName;
    for(Node* n : fixedTree.getDownPassSequence()){
        if(n->getIsTip()){
            for(const std::string& s : datRN){
                if(s == n->getName()){
                    nodeIdxByName.insert({s, n->getIndex()});
                    break;
                }
            }
        }
    }

    //----------Data wrangling-------//
    std::map<std::string, Eigen::MatrixXd> tipMatrices;
    for(int i = 0; i < datRN.size(); i++){
        std::string taxName = datRN[i];
        Eigen::MatrixXd datRow = dat->row(i);
        if (tipMatrices.find(taxName) == tipMatrices.end()) {
            tipMatrices.insert({taxName, datRow});
        } else {
            Eigen::MatrixXd& tipDat = tipMatrices[taxName];
            tipDat.conservativeResize(tipDat.rows() + 1, Eigen::NoChange);
            tipDat.row(tipDat.rows() - 1) = datRow;
        }
    }

    if(settings.getLogTransformData() == true)
        for(auto& s : tipMatrices)
            s.second = s.second.array().log();

    tipNames.reserve(tipMatrices.size());
    tipModels.reserve(tipMatrices.size());

    for(auto& s : tipMatrices){
        tipNames.push_back(s.first);

        TipModelV2* newTipModel = nullptr;
        if(settings.getWithIntraspecific() == false){
            Eigen::MatrixXd mean = Eigen::MatrixXd::Zero(1, s.second.cols());
            Eigen::VectorXi nObserved = Eigen::VectorXi::Zero(s.second.cols());
            for(int i = 0; i < s.second.rows(); i++){
                for(int j = 0; j < s.second.cols(); j++){
                    if(!std::isnan(s.second(i,j))){
                        mean(0, j) += s.second(i,j);
                        nObserved(j)++;
                    }
                }
            }
            for(int j = 0; j < s.second.cols(); j++){
                if(nObserved(j) == 0)
                    Msg::error("Taxon " + s.first + " has no observed values for trait " + std::to_string(j) + "; cannot compute a species mean");
                mean(0, j) /= (double)nObserved(j);
            }
            newTipModel = new TipModelV2(s.first, mean, this);
        }else{
            newTipModel = new TipModelV2(s.first, s.second, this);
        }
        tipModels.push_back(newTipModel);   // tipModels[k] now guaranteed to correspond to tipNames[k]
    }

    // resolve node indices in the same order as tipNames/tipModels
    tipNodeIdxs.resize(tipNames.size());
    for(int idx = 0; idx < (int)tipNames.size(); idx++){
        auto it = nodeIdxByName.find(tipNames[idx]);
        if(it == nodeIdxByName.end())
            Msg::error("Could not find tree node index for taxon " + tipNames[idx]);
        tipNodeIdxs[idx] = it->second;
    }

    tipMeansConcat.resize(tipModels.size(), dat->cols());
    for(int idx = 0; idx < (int)tipModels.size(); idx++)
        tipMeansConcat.row(idx) = tipModels[idx]->getTipMean();

    //mvBM set up
    MultivariateBrownianMotionV2::addData(tipNames, &tipMeansConcat);
    MultivariateBrownianMotionV2::setTree(&fixedTree);

    if((int)tipModels.size() != fixedTree.getNumTaxa())
        Msg::error("diff # of tip models instantiated from taxa");
}

PerikymataHSPv4::~PerikymataHSPv4(void){
    for(TipModelV2* t : tipModels)
        delete t;
}

double PerikymataHSPv4::lnLikelihood(void){
    //CANNOT CACHE LNL BECAUSE TIP MEANS CHANGE
    double lnl = 0.0;

    if(withPhylo== true)
        lnl += MultivariateBrownianMotionV2::lnLikelihood();

    for(TipModelV2* t : tipModels)
        lnl += t->lnLikelihood();
    return lnl;
}


double PerikymataHSPv4::lnPriorProbability(void){
    double lnp = MultivariateBrownianMotionV2::lnPriorProbability();
    for(TipModelV2* t : tipModels)
        lnp += t->lnPriorProbability();
    return lnp;
}

std::vector<std::string> PerikymataHSPv4::getParameterNames(void){
    std::vector<std::string> parmValues = MultivariateBrownianMotionV2::getParameterNames();
    for(TipModelV2* t : tipModels){
        std::vector<std::string> n = t->getParameterNames();
        parmValues.insert(parmValues.end(), n.begin(), n.end());
    }
    return parmValues;
}

std::vector<double> PerikymataHSPv4::getParameterString(void){
    parmValues = MultivariateBrownianMotionV2::getParameterString();
    for(TipModelV2* t : tipModels){
        scratchVec = t->getParameterString();
        parmValues.insert(parmValues.end(), scratchVec.begin(), scratchVec.end());
    }
    return parmValues;
}

void PerikymataHSPv4::print(void){
    for(TipModelV2* t : tipModels)
        t->print();
}

double PerikymataHSPv4::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    if(rng.uniformRv() < 0.9 && updateTipsOn == true){
        tipUpdate = true;
        
        int drawIdx = (int)(rng.uniformRv() * tipModels.size());
        if(drawIdx >= (int)tipModels.size())
            drawIdx = (int)tipModels.size() - 1;

        updatedTipModel = tipModels[drawIdx];
        savedTipIdx     = tipNodeIdxs[drawIdx];
        savedTipDataRow = MultivariateBrownianMotionV2::tipData.row(savedTipIdx);

        const double hr = updatedTipModel->update();
        if(hr == std::numeric_limits<double>::max())
            return hr;

        MultivariateBrownianMotionV2::tipData.row(savedTipIdx) = updatedTipModel->getTipMean();
        MultivariateBrownianMotionV2::instantiateIndependentContrasts();
        return hr;
    }else{
        tipUpdate = false;
        double hr = MultivariateBrownianMotionV2::update();
        return hr;
    }
}

void PerikymataHSPv4::updateForAcceptance(void){
    if(tipUpdate == false)
        MultivariateBrownianMotionV2::updateForAcceptance();
    else
        updatedTipModel->updateForAcceptance();
}

void PerikymataHSPv4::updateForRejection(void){
    if (tipUpdate == false) {
        MultivariateBrownianMotionV2::updateForRejection();
    } else {
        updatedTipModel->updateForRejection();
        MultivariateBrownianMotionV2::tipData.row(savedTipIdx) = savedTipDataRow;
        MultivariateBrownianMotionV2::instantiateIndependentContrasts();
    }
}
