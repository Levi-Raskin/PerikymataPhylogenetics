#include "Node.hpp"
#include "Msg.hpp"
#include "MultivariateBrownianMotionV2.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterTree.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#include <iostream>

MultivariateBrownianMotionV2::MultivariateBrownianMotionV2(void) : PhylogeneticModel(),
    fixedTree(false),
    branchLengthsInstantiated(false),
    log2Pi(std::log(2.0 * M_PI)),
    cachedLnL(-std::numeric_limits<double>::infinity()),
    cachedLnP(-std::numeric_limits<double>::infinity()){
    
}

void MultivariateBrownianMotionV2::addData(std::vector<std::string> rn, Eigen::MatrixXd* data){
    originalData = data;
    originalDataRownames = rn;
        
    // ─`─ Data dimensions ──────────────────────────────────────────────────────
    numberOfTraits = (int)originalData->cols();
    numberOfTips   = (int)originalData->rows();

    // ── Tree setup ───────────────────────────────────────────────────────────
    treeParam = new ParameterTree(1.0, originalDataRownames, 10.0);
    treeParam->forceBinary();
    parameters.push_back(treeParam);
    treeParam->getTree()->initializeDownPassSequence();
    dpseq             = treeParam->getTree()->getDownPassSequence();
    numberOfNodes     = (int)dpseq.size();
    numberOfInternalNodes = numberOfNodes - numberOfTips;

    // ── Node-indexed storage ─────────────────────────────────────────────────
    int maxNodeIndex = 0;
    for(Node* n : dpseq)
        maxNodeIndex = std::max(maxNodeIndex, n->getIndex());
    branchLength.resize(maxNodeIndex + 1);
    nodeVals.resize(maxNodeIndex + 1);
    nDesc.resize(2);
    modifiedBranches.reserve(dpseq.size());

    // ── Tip data ─────────────────────────────────────────────────────────────
    tipData.resize(numberOfTips, numberOfTraits);
    populateObservedData();

    // ── Contrast computation scratch space ───────────────────────────────────
    contrastScatterMatrix.resize(numberOfTraits, numberOfTraits);
    contrasts.reserve(numberOfInternalNodes);

    // ── IW prior / posterior ─────────────────────────────────────────────────
    dof = numberOfTraits + 2; //such that the mean is the scale matrix
    psi  = Eigen::MatrixXd::Identity(numberOfTraits, numberOfTraits);
    varianceCovarianceMatrix = Eigen::MatrixXd::Identity(numberOfTraits, numberOfTraits);
    mu0 = Eigen::VectorXd::Zero(numberOfTraits);
    dofN = dof + numberOfInternalNodes;

    // ── Finalize ─────────────────────────────────────────────────────────────
    normalizeProposalProbability();
    updateVarianceCovarianceMatrix();
}


MultivariateBrownianMotionV2::~MultivariateBrownianMotionV2(void){
    delete treeParam;
}

double MultivariateBrownianMotionV2::calculatePosteriorProbability(void){
    return Probability::InverseWishart::lnPdf(varianceCovarianceMatrix, psiN, dofN);
}

std::vector<std::string> MultivariateBrownianMotionV2::getParameterNames(void){
    std::vector<std::string> parmValues;
    
    //variance covariance matrix
    for(int i = 0; i < numberOfTraits; i++)
        for(int j = 0; j < numberOfTraits; j++)
            parmValues.push_back("vcv" + std::to_string(i) + "," +  std::to_string(j));

    return parmValues;
}

std::vector<double> MultivariateBrownianMotionV2::getParameterString(void){
    std::vector<double> parmValues;
    for(int i = 0; i < numberOfTraits; i++)
        for(int j = 0; j < numberOfTraits; j++)
            parmValues.push_back(varianceCovarianceMatrix(i,j));
    return parmValues;
}

double MultivariateBrownianMotionV2::lnLikelihood(void){
//    Eigen::MatrixXd SigmaInv = varianceCovarianceMatrix.inverse();
//    double logDetSigma = std::log(varianceCovarianceMatrix.determinant());
//    
//    cachedLnL = -0.5 * numberOfInternalNodes * (numberOfTraits * log2Pi + logDetSigma);
//    for(const Eigen::VectorXd& c : contrasts)
//        cachedLnL -= 0.5 * c.dot(SigmaInv * c);
//
//    return cachedLnL;
    Eigen::LLT<Eigen::MatrixXd> llt(varianceCovarianceMatrix);
    const Eigen::MatrixXd& L = llt.matrixL();

    double logDetSigma = 2.0 * L.diagonal().array().log().sum();

    Eigen::MatrixXd C(numberOfTraits, contrasts.size());
    for(int i = 0; i < (int)contrasts.size(); i++)
        C.col(i) = contrasts[i];

    Eigen::MatrixXd Z = llt.matrixL().solve(C);
    double traceSum = Z.squaredNorm();

    cachedLnL = -0.5 * numberOfInternalNodes * (numberOfTraits * log2Pi + logDetSigma) - 0.5 * traceSum;
    return cachedLnL;
}

void MultivariateBrownianMotionV2::instantiateIndependentContrasts(void){
    if(updatedParameter == treeParam || branchLengthsInstantiated == false ){
        Tree* t = treeParam->getTree();
        t->initializeDownPassSequence();
        dpseq = t->getDownPassSequence();
        root = t->getRoot();
                
        // Populate branch lengths
        for(Node* n : dpseq) {
            if(n != root) {
                Node* nAnc = n->getAncestor();
                branchLength[n->getIndex()] = t->getBranchLength(n, nAnc);
            }
        }
        
        branchLengthsInstantiated = true;
    }

    modifiedBranches.clear();
    contrasts.clear();
    
    //Instantiate independent contrasts
    
    for(Node* n : dpseq){
        int nIdx = n->getIndex();
        
        if(n->getIsTip() == false){
            nDesc = n->getDescendants();
            const Eigen::VectorXd& x0 = nodeVals[nDesc[0]->getIndex()];
            const Eigen::VectorXd& x1 = nodeVals[nDesc[1]->getIndex()];
            double v0 = branchLength[nDesc[0]->getIndex()];
            double v1 = branchLength[nDesc[1]->getIndex()];
            double blSum = v0 + v1;
            double invBlSum = 1.0 / blSum;
                
            contrasts.push_back((x1 - x0) / std::sqrt(blSum)); //felsensteins standardized contrast
            
            nodeVals[nIdx] = ( v1*x0 + v0*x1 ) * invBlSum; //this is the pruned node estimate;
            
            if(n != root){
                double origBL = branchLength[nIdx];
                modifiedBranches.emplace_back(nIdx, origBL);
                branchLength[nIdx] = origBL + (v0 * v1) * invBlSum;
            }
            
        }else{
            nodeVals[nIdx] = tipData.row(nIdx);
        }
    }
    
    for(auto& [idx, origVal] : modifiedBranches)
        branchLength[idx] = origVal;
}

double MultivariateBrownianMotionV2::lnPriorProbability(void){
    cachedLnP = treeParam->lnProbability();
    cachedLnP += Probability::InverseWishart::lnPdf(&varianceCovarianceMatrix, &psi, dof);
    return cachedLnP;
}

void MultivariateBrownianMotionV2::normalizeProposalProbability(void){
    double sum = 0.0;
    for (Parameter* p : parameters)
        sum += p->getProposalProbability();
    if(sum > 0)
        for (Parameter* p : parameters)
            p->setProposalProbability(p->getProposalProbability()/sum);
}

void MultivariateBrownianMotionV2::populateObservedData(void){
    
    for(int j = 0; j < originalDataRownames.size(); j++){
        int scratchIdx = -1;
        //find the tip associated with s
        const std::string& s = originalDataRownames[j];
        for(int i = 0; i < numberOfNodes; i++){
            Node* n = dpseq[i];
            if(n->getName() == s){
                scratchIdx = n->getIndex();
                tipData.row(scratchIdx) = originalData->row(j);
                break;
            }
        }
        if(scratchIdx == -1)
            Msg::error("Could not find tip " + s + " in tree");
    }
    
//    if(numberOfTips != treeParam->getTree()->getNumTaxa())
//        Msg::error("Data has disequal number of rows as tips");
}

void MultivariateBrownianMotionV2::print(void){
    for(Parameter* p : parameters)
        if(p->getParmPrintConsole() == true)
            std::cout << p->getName() << ": " << p->getAcceptanceRatio() << " " << p->getAdaptiveProposalActive() << " | ";
    std::cout << std::endl;
}

void MultivariateBrownianMotionV2::setTree(Tree* t){
    treeParam->setTree(t);
    treeParam->setProposalProbability(0.0);
    treeParam->setParmPrintConsole(false);
    
    normalizeProposalProbability();
            
    fixedTree = true;
    treeParam->getTree()->initializeDownPassSequence();
    dpseq = t->getDownPassSequence();
    
    int maxNodeIndex = 0;
    for(Node* n : dpseq)
        maxNodeIndex = std::max(maxNodeIndex, n->getIndex());
    branchLength.resize(maxNodeIndex + 1);
    nodeVals.resize(maxNodeIndex + 1);
    
    populateObservedData(); //need to repopulate because of edgecase where tip idxs changed
    updateVarianceCovarianceMatrix();
}

double MultivariateBrownianMotionV2::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hr  = 0.0;
    if(fixedTree == false && rng.uniformRv() < 0.5){
        updatedParameter = treeParam; //update tree
        hr = updatedParameter->update();
    }else{
        updatedParameter = nullptr; // variance covariance matrix
        updateVarianceCovarianceMatrix();
        hr = std::numeric_limits<double>::max();
    }
    instantiateIndependentContrasts();
    return hr;
}

void MultivariateBrownianMotionV2::updateForAcceptance(void){
    if(updatedParameter != nullptr)
        updatedParameter->updateForAcceptance();
}

void MultivariateBrownianMotionV2::updateForRejection(void){
    if(updatedParameter != nullptr)
        updatedParameter->updateForRejection();
}

void MultivariateBrownianMotionV2::updateVarianceCovarianceMatrix(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    contrastScatterMatrix = Eigen::MatrixXd::Zero(numberOfTraits, numberOfTraits);
    for(const Eigen::VectorXd& c : contrasts)
        contrastScatterMatrix.noalias() += c * c.transpose();
    
    psiN = psi + contrastScatterMatrix;
    Eigen::MatrixXd psiNInvLower = psiN.inverse().llt().matrixL();
    
    Probability::InverseWishart::rv(&rng, varianceCovarianceMatrix, psiNInvLower, dofN);
}
