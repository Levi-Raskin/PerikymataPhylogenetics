#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "ReadTSV.hpp"
#include "SimulateData.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "Utility.hpp"

#include <iostream>
#include <fstream>

SimulateData::SimulateData(void) : tree(nullptr), trials(0){
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
    trials++;
    rownames.clear();
    trueTipMeans.clear();
    trueTipVCVs.clear();
    tipNameToIndex.clear();
    trueMissingValues.clear();

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    //---simulate tree---//
    if(tree != nullptr)
        delete tree;
    tree = new Tree(tipNames, 10.0);
//    tree = new Tree("(((Gorilla_beringei:2.558516,Gorilla_gorilla:2.558516):6.093717,((Homo_sapiens:0.568721,Neanderthal:0.538721):5.607159,(Pan_paniscus:2.333553,Pan_troglodytes:2.333553):3.842326):2.476353):6.480222,(Pongo_abelii:3.825854,Pongo_pygmaeus:3.825854):11.306601);");
//    ntips = tree->getNumTaxa();
        
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
        Eigen::MatrixXd tipVCV = Probability::InverseWishart::rv(&rng, psi, priorDOF);

        std::string tipName = taxNames[i];
        Eigen::VectorXd tipMean = dat.row(i);
        
        std::vector<std::string> colnames = {};
        for(int j =0; j < ntraits; j++)
            colnames.push_back("trait" + std::to_string(j));

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

void SimulateData::checkCredInt(void){
    UserSettings& settings = UserSettings::userSettings();
    ReadTSV r(settings.getOutputFile() + "Outfile.tsv", false, true);
    Eigen::MatrixXd rDat = r.getEigenMat();
    std::vector<std::string> cn = r.getColnames();
    
    // Evaluate evolutionary VCV coverage
    int evoVCVCovered = 0;
    for(int x = 0; x < cn.size(); x++){
        std::string s = cn[x];
        if(s.substr(0, 7) == "evo_vcv"){
            size_t rowPos = s.find_first_of("0123456789");
            size_t commaPos = s.find(',');
            int i = std::stoi(s.substr(rowPos, commaPos - rowPos));
            int j = std::stoi(s.substr(commaPos + 1));
            auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.1);
            if(sampledEvoVCV(i,j) < impInterval.second && sampledEvoVCV(i,j) > impInterval.first){
                vcvInCredInt(i,j)++;
                evoVCVCovered++;
            }
        }
    }
    // Evaluate tip mean coverage
    for(auto& tipEntry : trueTipMeans) {
        std::string tipName = tipEntry.first;
        Eigen::VectorXd trueMean = tipEntry.second;
        int tipIdx = tipNameToIndex[tipName];
        int tipMeanCovered = 0;
        
        for(int traitIdx = 0; traitIdx < ntraits; traitIdx++) {
            std::string paramName = tipName + "_mean_" + std::to_string(traitIdx);
            
            auto it = std::find(cn.begin(), cn.end(), paramName);
            if(it != cn.end()) {
                int colIdx = std::distance(cn.begin(), it);
                auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(colIdx), 0.1);
                
                if(trueMean(traitIdx) < impInterval.second && trueMean(traitIdx) > impInterval.first) {
                    tipMeanInCredInt(tipIdx,traitIdx)++;
                    tipMeanCovered++;
                }
            }
        }
    }
    
    // Evaluate tip VCV coverage
    for(auto& tipEntry : trueTipVCVs) {
        std::string tipName = tipEntry.first;
        Eigen::MatrixXd trueVCV = tipEntry.second;
        int tipIdx = tipNameToIndex[tipName];
        int tipVCVCovered = 0;
        
        for(int i = 0; i < ntraits; i++) {
            for(int j = 0; j < ntraits; j++) {
                std::string paramName = tipName + "_vcv_(" + std::to_string(i) + "," + std::to_string(j) + ")";
                
                auto it = std::find(cn.begin(), cn.end(), paramName);
                if(it != cn.end()) {
                    int colIdx = std::distance(cn.begin(), it);
                    auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(colIdx), 0.1);
                    
                    if(trueVCV(i,j) < impInterval.second && trueVCV(i,j) > impInterval.first) {
                        tipVCVInCredInt[tipIdx](i,j)++;
                        tipVCVCovered++;
                    }
                }
            }
        }
    }
    
    // Evaluate missing data imputation
    for(int x = 0; x < (int)cn.size(); x++){
        std::string s = cn[x];
        if(s.substr(0, 8) == "missing_"){
            size_t parenOpen  = s.find('(');
            size_t parenClose = s.find(')');
            size_t lastUnderscore = s.rfind('_', parenOpen - 1);

            std::string tipName = s.substr(8, lastUnderscore - 8);

            std::string inside = s.substr(parenOpen + 1, parenClose - parenOpen - 1);
            size_t commaPos = inside.find(',');
            int localRow = std::stoi(inside.substr(0, commaPos));
            int col      = std::stoi(inside.substr(commaPos + 1));

            // reconstruct absolute row
            auto tipIt = tipNameToIndex.find(tipName);
            if(tipIt == tipNameToIndex.end()) continue;
            int tipIdx = tipIt->second;
            int absRow = tipIdx * nind + localRow;

            auto it = trueMissingValues.find({absRow, col});
            if(it != trueMissingValues.end()){
                double trueVal = it->second;
                auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.1);
                if(trueVal > impInterval.first && trueVal < impInterval.second){
                    int impIdx = (int)std::distance(trueMissingValues.begin(), it);
                    imputedInCredInt(impIdx)++;
                }
            }
        }
    }
}

void SimulateData::print(void){
    int total = 0;
    for(auto& m : tipVCVInCredInt)
        total += m.sum();

    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Evolutionary VCV coverage:            " << vcvInCredInt.sum() << "/" << (trials * ntraits * ntraits) << "\t\t | (" << (double)vcvInCredInt.sum() / (trials * ntraits * ntraits) << ")" << "\n";
    std::cout << "Tip VCV coverage:                     " << total << "/" << (trials * ntips * ntraits * ntraits) << "\t\t | (" << (double)total / (trials * ntips * ntraits * ntraits) << ")" << "\n";
    std::cout << "Tip mean coverage:                    " << tipMeanInCredInt.sum() << "/" << (trials * ntips * ntraits) << "\t\t | (" << (double)tipMeanInCredInt.sum() / (trials * ntips * ntraits) << ")" << "\n";
    std::cout << "Missing data coverage:                " << imputedInCredInt.sum() << "/" << (trials * nimp) << "\t\t | (" << (double)imputedInCredInt.sum() / (trials * nimp) << ")" << "\n";
    std::cout << "-----------------------------------------------------------------------" << std::endl;
}

void SimulateData::writeCoverage(void){
    UserSettings& settings = UserSettings::userSettings();
    std::string logFile = settings.getOutputFile() + "CoverageResults.txt";
    
    std::ofstream log(logFile);
    if (!log.is_open())
        Msg::error("Could not open log file: " + logFile);
        
    int total = 0;
    for(auto& m : tipVCVInCredInt)
        total += m.sum();
    
    if(settings.getWithIntraspecific() == false){
        log << "Total coverage: " << vcvInCredInt.sum() << "/" << (trials * ntraits * ntraits) << " | (" << (double)vcvInCredInt.sum() / (trials * ntraits * ntraits) << ")\n";
    }else if(settings.getWithPhylogeny() == false){
        int cumCov = total;
        cumCov += tipMeanInCredInt.sum();
        cumCov += imputedInCredInt.sum();
        log << "Total coverage: " << cumCov << "/" << ((trials * ntips * ntraits * ntraits) + (trials * ntips * ntraits) + (trials * nimp)) << " | ("<< (double) cumCov / ((trials * ntips * ntraits * ntraits) + (trials * ntips * ntraits) + (trials * nimp)) << ")\n";
    }else{
        int cumCov = total;
        cumCov += vcvInCredInt.sum();
        cumCov += tipMeanInCredInt.sum();
        cumCov += imputedInCredInt.sum();

        log << "Total coverage: " << cumCov << "/" << ((trials * ntraits * ntraits) + (trials * ntips * ntraits * ntraits) + (trials * ntips * ntraits) + (trials * nimp)) << " | ("<< (double) cumCov / ((trials * ntraits * ntraits) + (trials * ntips * ntraits * ntraits) + (trials * ntips * ntraits) + (trials * nimp)) << ")\n";
    }
    log << "-----------------------------------------------------------------------\n";
    log << "Evolutionary VCV coverage:            " << vcvInCredInt.sum() << "/" << (trials * ntraits * ntraits) << " | (" << (double)vcvInCredInt.sum() / (trials * ntraits * ntraits) << ")" << "\n";
    log << "Tip VCV coverage:                     " << total << "/" << (trials * ntips * ntraits * ntraits) << " | (" << (double)total / (trials * ntips * ntraits * ntraits) << ")" << "\n";
    log << "Tip mean coverage:                    " <<  tipMeanInCredInt.sum() << "/" << (trials * ntips * ntraits) << " | (" << (double)tipMeanInCredInt.sum() / (trials * ntips * ntraits) << ")" << "\n";
    log << "Missing data coverage:                " <<  imputedInCredInt.sum() << "/" << (trials * nimp) << " | (" << (double)imputedInCredInt.sum() / (trials * nimp) << ")" << "\n";

    log.close();

}
