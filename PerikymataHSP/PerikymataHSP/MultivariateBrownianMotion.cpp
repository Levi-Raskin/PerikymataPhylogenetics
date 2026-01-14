#include "Node.hpp"
#include "Msg.hpp"
#include "MultivariateBrownianMotion.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterTree.hpp"
#include "ParameterInvVCV.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TicToc.hpp"
#include "Utility.hpp"

#include <iostream>
#include <random>

MultivariateBrownianMotion::MultivariateBrownianMotion(std::vector<std::string> rn, Eigen::MatrixXd* data) : PhylogeneticModel(){
    
    addData(rn, data);
    
    numberOfTraits = (int)observedData->cols();
    varianceCovarianceMatrix = new ParameterInvVCV(1.0, "vcvInv", numberOfTraits, this);
    parameters.push_back( varianceCovarianceMatrix );
    
    ParameterTree* pt = new ParameterTree(1.0, this, observedDataRownames, 10.0);
    parameters.push_back( pt );
    
    double sum = 0.0;
    for (Parameter* p : parameters)
        sum += p->getProposalProbability();
    for (Parameter* p : parameters)
        p->setProposalProbability(p->getProposalProbability()/sum);
}

void MultivariateBrownianMotion::addData(std::vector<std::string> rn, Eigen::MatrixXd* data){
    observedData = data;
    Eigen::MatrixXd X_centered = *observedData;
    Eigen::RowVectorXd empiricalMean = X_centered.colwise().mean();
    X_centered.rowwise() -= empiricalMean;
    datScatter = X_centered.transpose() * X_centered;
    
    observedDataRownames = rn;
    tipMeanOrdering = rn;
    
    //NEW//
    for (int i = 0; i < tipMeanOrdering.size(); ++i)
        tipHashMap[tipMeanOrdering[i]] = i;
}

std::vector<std::string> MultivariateBrownianMotion::getParameterNames(void){
    std::vector<std::string> parmValues;
    for(int p = 0; p < parameters.size(); p++){
        ParameterMatrix* pt = dynamic_cast<ParameterMatrix*>(parameters[p]);
        if(pt != nullptr){
            Eigen::MatrixXd scratch = (pt->getValue());
            for(int i = 0; i < scratch.rows(); i++)
                for(int j = 0; j < scratch.cols(); j++)
                    parmValues.push_back(pt->getName() + std::to_string(i) + "," +  std::to_string(j));
        }else
            parmValues.push_back(parameters[p]->getName());
    }
    parmValues.push_back("adaptiveProposalTuningAcive");
    return parmValues;

}

std::vector<double> MultivariateBrownianMotion::getParameterString(void){
    std::vector<double> parmValues;
    for(int i = 0; i < parameters.size(); i++){
        ParameterDouble* pt = dynamic_cast<ParameterDouble*>(parameters[i]);
        ParameterMatrix* mat = dynamic_cast<ParameterMatrix*>(parameters[i]);
        
        if(mat != nullptr){
            Eigen::MatrixXd scratch = (mat->getValue());
            for(int i = 0; i < scratch.rows(); i++)
                for(int j = 0; j < scratch.cols(); j++)
                    parmValues.push_back(scratch(i,j));
        }else if(pt != nullptr){
            parmValues.push_back(pt->getValue());
        }else{
            parmValues.push_back(-1.0);
        }
    }
    bool aProp = false;
    for(Parameter* p : parameters)
        if(p->getAdaptiveProposalActive() == true){
            aProp = true;
            break;
        }
    parmValues.push_back(aProp);
    return parmValues;
}

double MultivariateBrownianMotion::lnLikelihood(void){
    double lnL = lnLikelihoodF85Optim();
    return lnL;
}

double MultivariateBrownianMotion::lnLikelihoodF85Optim(void){
    //Expects varianceCovarianceMatrix to be inverse
//    Msg::warning("Eigen solver implementation is active");
    Tree pruningTree = *(getTree());
    pruningTree.initializeDownPassSequence();
    std::vector<Node*> dpseq = pruningTree.getDownPassSequence();

    double lnL = 0.0;
    
    //Precalcualtions
    const double log2Pi = std::log(2 * M_PI);
    Eigen::MatrixXd vcvInv = varianceCovarianceMatrix->getValue();
    double halfLogDet = 0.5 * (std::log(1) - std::log(vcvInv.determinant()));
//    if(std::isnan(halfLogDet)){
//        Utility::EigenUtils::printEigenR(vcvInv.inverse());
//        Msg::error("here");
//    }
    const double halfNTraits = -numberOfTraits/2;
    const double halfNTraitsLog2PiProd = halfNTraits * log2Pi;
    
    //Preallocations
    std::unordered_map<Node*, Eigen::VectorXd> nodeVals;
    std::vector<Node*> nDesc;
    nDesc.resize(2);
    Eigen::VectorXd u1 = Eigen::VectorXd::Zero(numberOfTraits);
    Eigen::VectorXd extantMeans = Eigen::VectorXd::Zero(numberOfTraits);
    
    for(Node* n : dpseq){
        if(n->getIsTip() == false){
            nDesc = n->getDescendants();
            const Eigen::VectorXd& x0 = nodeVals[nDesc[0]];
            const Eigen::VectorXd& x1 = nodeVals[nDesc[1]];
            double v0 = nDesc[0]->getBranchLength();
            double v1 = nDesc[1]->getBranchLength();
            double blSum = v0 + v1;
            u1 = x1 - x0; //this is the raw contrast
            
            nodeVals[n] = ( v1*x0 + v0*x1 ) / (blSum); //this is the pruned node estimate;
            if(n != pruningTree.getRoot()){
                double currBL = n->getBranchLength();
                currBL += (v0 * v1) / (blSum);
                n->setBranchLength(currBL);
            }
            
            lnL += halfNTraitsLog2PiProd + halfNTraits * std::log(blSum) + halfLogDet - (1/(2*(blSum))) * u1.dot(vcvInv * u1);
        }else
            nodeVals[n] = (*observedData).row(tipHashMap.at(n->getName()));
    }
    return lnL;
    
    //The following matches the likelihood calc in Alvarez-Carretero 2019
    //Confirmed identical to the old, working, unoptim code
//    Tree pruningTree = *(getTree());
//    pruningTree.initializeDownPassSequence();
//    std::vector<Node*>& dpseq = pruningTree.getDownPassSequence();
//
//    double lnL = 0.0;
//
//    //Precalcualtions
//    double log2Pi = std::log(2 * M_PI);
//    Eigen::MatrixXd vcv = varianceCovarianceMatrix->getValue();
//    double det = vcv.determinant();
//    Eigen::MatrixXd vcvInv = vcv.inverse();
//    double halfLogDet = -0.5*std::log(det);
//    double halfNTraits = -numberOfTraits/2;
//
//    //Preallocations
//    std::unordered_map<Node*, Eigen::VectorXd> nodeVals;
//    std::vector<Node*> nDesc;
//    nDesc.resize(2);
//    Eigen::VectorXd u1 = Eigen::VectorXd::Zero(numberOfTraits);
//    Eigen::VectorXd extantMeans = Eigen::VectorXd::Zero(numberOfTraits);
//
//    for(Node* n : dpseq){
//        if(n->getIsTip() == false){
//            nDesc = n->getDescendants();
//            const Eigen::VectorXd& x0 = nodeVals[nDesc[0]];
//            const Eigen::VectorXd& x1 = nodeVals[nDesc[1]];
//            double v0 = pruningTree.getBranchLength(n, nDesc[0]);
//            double v1 = pruningTree.getBranchLength(n, nDesc[1]);
//            double blSum = v0 + v1;
//            u1 = x1 - x0; //this is the raw contrast
//
//            nodeVals[n] = ( v1*x0 + v0*x1 ) / (blSum); //this is the pruned node estimate;
//
//            if(n != pruningTree.getRoot()){
//                double currBL = pruningTree.getBranchLength(n, n->getAncestor());
//                currBL += (v0 * v1) / (blSum);
//                pruningTree.setBranch(n, n->getAncestor(), currBL);
//            }
//
//            lnL += halfNTraits * (log2Pi +  std::log(blSum)) + halfLogDet - (1/(2*(blSum))) * (u1.transpose() * (vcvInv * u1)).value();
//        }else
//            nodeVals[n] = (*observedData).row(tipHashMap.at(n->getName()));
//    }
//    return lnL;
    
    //LYR old, unoptim, but WORKING code
    /*
     Tree pruningTree = *(getTree());
     pruningTree.checkBranchLengthsNeg();
     pruningTree.initializeDownPassSequence();
     std::vector<Node*> dpseq = pruningTree.getDownPassSequence();

     double lnL = 0.0;
     
     double log2Pi = std::log(2 * M_PI);
     std::map<Node*, Eigen::VectorXd> nodeVals;
     for(Node* n : dpseq){
         if(n->getIsTip() == false){
             std::vector<Node*> nDesc = n->getDescendants();
             if(nDesc.size() != 2){
                 pruningTree.print();
                 std::cout << "Issue node: " << n->getIndex() << std::endl;
                 Msg::error("Expecting 2 decendants for all interior nodes");
             }
             Eigen::VectorXd x0 = nodeVals[nDesc[0]];
             Eigen::VectorXd x1 = nodeVals[nDesc[1]];
             double v0 = pruningTree.getBranchLength(n, nDesc[0]);
             double v1 = pruningTree.getBranchLength(n, nDesc[1]);
             
             Eigen::VectorXd u1 = x1 - x0; //this is the raw contrast

             Eigen::VectorXd xInt = ( v1*x0 + v0*x1 ) / (v0 + v1); //this is the pruned node estimate;
             nodeVals.insert({n, xInt});
             
             //now we need to lengthen the branch below n
             //this is lengthed by the variance of product of two normals
             if(n != pruningTree.getRoot()){
                 double currBL = pruningTree.getBranchLength(n, n->getAncestor());
                 currBL += (v0 * v1) / (v0 + v1);
                 pruningTree.setBranch(n, n->getAncestor(), currBL);
             }
 //            Eigen::VectorXd mean = Eigen::VectorXd::Zero(u1.size());
 //            Eigen::MatrixXd vcv = (varianceCovarianceMatrix->getValue())*(v0 + v1);
             
             double det = (varianceCovarianceMatrix->getValue()).determinant();
             Eigen::MatrixXd vcvInv = varianceCovarianceMatrix->getValue().inverse();

             lnL += -numberOfTraits/2 * (log2Pi +  std::log(v0 + v1)) + -0.5*std::log(det) + (-(1/(2*(v0+v1))) * u1.transpose() * vcvInv * u1);
             if(std::isnan(lnL))
                 Msg::error("lnL is nan");
         }else{
             int rowIdx = -1;
             for(int i = 0; i < tipMeanOrdering.size(); i++){
                 if(n->getName() == tipMeanOrdering[i]){
                     rowIdx = i;
                     break;
                 }
             }
             if(rowIdx == -1)
                 Msg::error("did not find tip");
             
             Eigen::VectorXd extantMeans = Eigen::VectorXd::Zero(numberOfTraits);
             for(int i = 0; i < numberOfTraits; i++)
                 extantMeans(i) = (*observedData)(rowIdx, i);
             
             nodeVals.insert({n, extantMeans});
         }
     }
     */
}

double MultivariateBrownianMotion::lnPriorProbability(void){
    double priorProb = 0.0;
    for(auto p : parameters)
        priorProb += p->lnProbability();
    return priorProb;
}

void MultivariateBrownianMotion::print(void){
    std::cout << " -- ";
    for(Parameter* p : parameters)
            std::cout << p->getName() << ": " << p->getAcceptanceRatio() << " " << p->getAdaptiveProposalActive() << " | ";
    std::cout << std::endl;
}

void MultivariateBrownianMotion::setTree(Tree* t){
    ParameterTree* parmTree = nullptr;
    for (Parameter* p : parameters){
        ParameterTree* pt = dynamic_cast<ParameterTree*>(p);
        if (pt != nullptr){
            parmTree = pt;
            break;
        }
    }
    if(parmTree != nullptr){
        parmTree->setTree(t);
        parmTree->setProposalProbability(0.0);
        double sum = 0.0;
        for (Parameter* p : parameters)
            sum += p->getProposalProbability();
        for (Parameter* p : parameters)
            p->setProposalProbability(p->getProposalProbability()/sum);
    }
}

double MultivariateBrownianMotion::update(void){
    Parameter* parm = nullptr;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
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
    return updatedParameter->update();
}

void MultivariateBrownianMotion::updateForAcceptance(void){
    updatedParameter->updateForAcceptance();
}

void MultivariateBrownianMotion::updateForRejection(void){
    updatedParameter->updateForRejection();
}
