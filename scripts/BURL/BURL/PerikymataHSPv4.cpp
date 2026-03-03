#include "Msg.hpp"
#include "MultivariateBrownianMotionV2.hpp"
#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterIntraspecificMean.hpp"
#include "ParameterTree.hpp"
#include "PerikymataHSPv4.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TipModelV2.hpp"
#include "Tree.hpp"

#include <iostream>
#include <string>
#include <vector>

PerikymataHSPv4::PerikymataHSPv4(Tree* backbone, std::vector<std::string> datRN, Eigen::MatrixXd* dat) : MultivariateBrownianMotionV2(),
    fixedTree(*backbone),
    updateTipsOn(true){
    
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
    
    for(Node* n : fixedTree.getDownPassSequence()){
        if(n->getIsTip()){
            bool found = false;
            for(const std::string& s : datRN){
                if(s == n->getName()){
                    found = true;
                    tipIdxs.insert({s, n->getIndex()});
                    break;
                }
            }
        }
    }


    //----------Data wrangling-------//

    // Constructing tipDataIncomplete objects
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
            
//    for(auto& s : tipMatrices){
//        Msg::warning("these data are log transformed");
//        s.second = s.second.array().log(); //log transforming pk/mm so that I don't have to worry about boudning
//    }
    
    for(auto&s : tipMatrices){
        tipNames.push_back(s.first);
        TipModelV2* newTipModel = new TipModelV2(s.first, s.second, this);
        tipModels.insert({s.first, newTipModel});
    }
    
    tipMeansConcat.resize(tipModels.size(), dat->cols());
    int idx = 0;
    for(auto&s : tipModels){
        tipMeansConcat.row(idx) = s.second->getTipMean();
        idx++;
    }
    
    //mvBM set up
    MultivariateBrownianMotionV2::addData(tipNames, &tipMeansConcat);
    MultivariateBrownianMotionV2::setTree(&fixedTree);
    
    if(tipModels.size() != fixedTree.getNumTaxa())
        Msg::error("diff # of tip models instantiated from taxa");
}

PerikymataHSPv4::~PerikymataHSPv4(void){
    for(auto&s : tipModels){
        delete s.second;
    }
}

double PerikymataHSPv4::lnLikelihood(void){
    //CANNOT CACHE LNL BECAUSE TIP MEANS CHANGE
    double lnl = MultivariateBrownianMotionV2::lnLikelihood();
    for (auto s : tipModels)
        lnl += s.second->lnLikelihood();
    return lnl;
}


double PerikymataHSPv4::lnPriorProbability(void){
    double lnp = MultivariateBrownianMotionV2::lnPriorProbability();
    for (auto s : tipModels)
        lnp += s.second->lnPriorProbability();
    return lnp;
}

std::vector<std::string> PerikymataHSPv4::getParameterNames(void){
    std::vector<std::string> parmValues = MultivariateBrownianMotionV2::getParameterNames();
    //adding in tipModel names
    for(auto& s : tipModels){
        std::vector<std::string> n = s.second->getParameterNames();
        parmValues.insert(parmValues.end(), n.begin(), n.end());
    }
    
    return parmValues;
}

std::vector<double> PerikymataHSPv4::getParameterString(void){
    parmValues = MultivariateBrownianMotionV2::getParameterString();
    
    //adding in tipModel values
    for(auto& s : tipModels){
        scratchVec = s.second->getParameterString();
        parmValues.insert(parmValues.end(), scratchVec.begin(), scratchVec.end());
    }
    
    return parmValues;
}

void PerikymataHSPv4::print(void){
    MultivariateBrownianMotionV2::print();
    for(auto& s : tipModels)
        s.second->print();
}

double PerikymataHSPv4::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    if(rng.uniformRv() < 0.9 && updateTipsOn == true){
        tipUpdate = true;
        retry:
        std::string s = tipNames[(int)(rng.uniformRv() * tipNames.size())];
        updatedTipModel = tipModels[s];
        if(updatedTipModel == nullptr)
            goto retry;
        double hr = updatedTipModel->update();
        MultivariateBrownianMotionV2::tipData.row(tipIdxs[s]) = updatedTipModel->getTipMean();
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
    if(tipUpdate == false)
        MultivariateBrownianMotionV2::updateForRejection();
    else
        updatedTipModel->updateForRejection();
}
