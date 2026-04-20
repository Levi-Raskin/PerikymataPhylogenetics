#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "SimulateData.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "Utility.hpp"

SimulateData::SimulateData(void) : tree(nullptr){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    UserSettings& settings = UserSettings::userSettings();
    ntraits = settings.getNumTraits();
    
    //specify prior
    priorDOF = ntraits + 2;
    psi = Eigen::MatrixXd::Constant(ntraits, ntraits, 1e-6);
    psi.diagonal().array() = 1.0;
    psi = Eigen::VectorXd::Ones(ntraits).asDiagonal() * psi * Eigen::VectorXd::Ones(ntraits).asDiagonal();
    
    //preallocs
    nreps = settings.getNumReps();
    ntips = settings.getNumTips();
    nimp = settings.getNumImputed();
    nind = settings.getNumObserved();
    
    vcvInCredInt = Eigen::MatrixXi::Zero(ntraits, ntraits);
    tipMeanInCredInt = Eigen::MatrixXi::Zero(ntips, ntraits);

    for(int i = 0; i < ntips; i++)
        tipVCVInCredInt.push_back(Eigen::MatrixXi::Zero(ntraits, ntraits));
    
    imputedInCredInt = Eigen::VectorXi::Zero(nimp);
    
    for(int i = 0; i < ntips; i++)
        tipNames.push_back("t" + std::to_string(i));
}

void SimulateData::simulateData(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    //---simulate tree---//
    if(tree != nullptr)
        delete tree;
        
    tree = new Tree(tipNames, 10.0);
        
    //---simulate evo VCV---//
    sampledEvoVCV = psi;
    Probability::InverseWishart::rv(&rng, sampledEvoVCV, psi, priorDOF);
    
    //---simulate data---//
    
    //mvBM
    std::vector<double> rootMean;
    for(int i = 0; i < ntraits; i++)
        rootMean.push_back(Probability::Normal::rv(&rng));
    std::map<Node*, std::vector<double>> simdat;
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    for (auto i = dpseq.rbegin(); i != dpseq.rend(); i++){
        Node* n = *i;
        if(n == tree->getRoot()){
            simdat.insert({n,rootMean});
        }else{
            std::vector<double> ancMean =(simdat[n->getAncestor()]);
            Eigen::MatrixXd var = sampledEvoVCV * tree->getBranchLength(n, n->getAncestor());
            std::vector<double> draw = Probability::MultivariateNormal::rv(&rng, ancMean, &var);
            simdat.insert({n,draw});
        }
    }

    //filter out means
    std::vector<std::vector<double>> taxDat;
    std::vector<std::string> taxNames;
    for(auto a : simdat){
        if(a.first->getIsTip() == true){
            taxNames.push_back(a.first->getName());
            taxDat.push_back(a.second);
        }
    }

    //add tip VCVs and simulate data at tips
    //vector of tip means
    Eigen::MatrixXd dat = Utility::EigenUtils::vectorMatrix2Eigen(taxDat);
        
    //sample individuals at each tip
    data = Eigen::MatrixXd::Zero(ntips*nind, ntraits);
    int idx = 0;
    for(int i = 0; i < ntips; i++){
        //sample tip VCV
        Eigen::VectorXd mu0 = Eigen::VectorXd::Zero(ntraits);
        Eigen::MatrixXd psi = Eigen::MatrixXd::Constant(ntraits, ntraits, 1e-6);
        psi.diagonal().array() = 1.0;
        psi = Eigen::VectorXd::Ones(ntraits).asDiagonal() * psi * Eigen::VectorXd::Ones(ntraits).asDiagonal();
        Eigen::MatrixXd tipVCV = Probability::InverseWishart::rv(&rng, psi, priorDOF);

        std::string tipName = taxNames[i];
        Eigen::VectorXd tipMean = dat.row(i);
        
        std::vector<std::string> colnames = {};
        for(int i =0; i < ntraits; i++)
            colnames.push_back("trait" + std::to_string(i));

        // Store true values
        trueTipMeans[tipName] = tipMean;
        trueTipVCVs[tipName] = tipVCV;
        tipNameToIndex[tipName] = i;
        
        for(int j = 0; j < nind; j++){
            Eigen::VectorXd tmprow = Probability::MultivariateNormal::rv(&rng, tipMean, &tipVCV);
            data.row(idx) = tmprow;
            rownames.push_back(tipName);
            idx++;
        }
    }

    //---add missing elements---//
    for(int i =0; i < nimp; i++){
        int missingI, missingJ;
        int cnt = 0;
        do {
            missingI = (int)(rng.uniformRv() * data.rows());
            missingJ = (int)(rng.uniformRv() * data.cols());
            cnt++;
            if(cnt > 100)
                Msg::error("nimp probably higher than # of data elemnts");
        } while (trueMissingValues.count({missingI, missingJ}) > 0);

        trueMissingValues.insert({{missingI, missingJ}, data(missingI, missingJ)});
        data(missingI, missingJ) = std::numeric_limits<double>::quiet_NaN();
    }
}


/*
        void                                    checkCredInt(void);
        double                                  getVCVInCredInt(void);
        double                                  getTipMeanInCredInt(void);
        double                                  getTipVCVInCredInt(void);
        double                                  getImputedInCredInt(void);
        void                                    simulateData(void);
        
    private:
//coverage instance vars
        Eigen::MatrixXi                         vcvInCredInt;       // ntraits × ntraits
        Eigen::MatrixXi                         tipMeanInCredInt;   // ntips   × ntraits
        std::vector<Eigen::MatrixXi>            tipVCVInCredInt;    // ntips   × ntraits × ntraits
        Eigen::VectorXi                         imputedInCredInt;   // nimp
        
        //prior parameters
        double                                  priorDOF;
        Eigen::MatrixXd                         psi;
        
        //simulated parameters
        Eigen::VectorXd                         `trueMissingValues`;
        
        //misc
        int                                     nimp;
        int                                     nreps;
        int                                     ntips;
        int                                     ntraits;*/
