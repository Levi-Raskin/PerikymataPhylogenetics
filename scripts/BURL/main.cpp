#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <random>
#include <regex>

#include "BrownianMotion.hpp"
#include "HomininDisparityAnalysis.hpp"
#include "LkjOnionSampler.hpp"
#include "Mcmc.hpp"
#include "MetropolisCoupledMcmc.hpp"
#include "MultivariateBrownianMotion.hpp"
#include "MultivariateBrownianMotionV2.hpp"
#include "Node.hpp"
#include "PerikymataHSPv3.hpp"
#include "PerikymataHSPv4.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "ReadCSV.hpp"
#include "ReadTSV.hpp"
#include "TicToc.hpp"
#include "TipModel.hpp"
#include "TipModelV2.hpp"
#include "Tree.hpp"
#include "Utility.hpp"
#include "UserSettings.hpp"
#include "WriteTSV.hpp"

int main(int argc, const char* argv[]) {

   
    UserSettings& settings = UserSettings::userSettings();
    settings.initializeSettings(argc, argv);
    settings.print();

    #if 0
    // Multivariate Brownian Motion v2
    int ntraits = 10;
    int in = 0;
    int vcvInCredInt[ntraits][ntraits];
    for(int i = 0; i < ntraits; i++)
        for(int j = 0; j < ntraits; j++)
            vcvInCredInt[i][j] = 0;
    int nreps = 100;
//    #pragma omp parallel for num_threads(10) shared(vcvInCredInt)
        for(int cyc = 0; cyc < nreps; cyc++){
            RandomVariable& rng = RandomVariable::randomVariableInstance();
            
            Eigen::MatrixXd psi = Eigen::MatrixXd::Identity(ntraits, ntraits);;
            Eigen::MatrixXd sampledVCV = psi;
            Probability::InverseWishart::rv(&rng, sampledVCV, psi, ntraits+2);

            #pragma omp critical
            {
                std::cout << "========================================" << "\n";
                std::cout << "Evolutionary variance covariance matrix for thread " << omp_get_thread_num() << "\n";
                Utility::EigenUtils::printEigen(sampledVCV);
                std::cout << "========================================" << "\n";
            }
            std::vector<std::string> taxnames;
            for(int i = 0; i < 50; i++)
                taxnames.push_back("t"+std::to_string(i));
            Tree t = Tree(taxnames, 10);
            t.forceBinary();
            t.scaleTreeHeight(1.0);
            
            //simulate data under mvBM
            std::vector<double> rootMean;
            for(int i = 0; i < sampledVCV.rows(); i++)
                rootMean.push_back(Probability::Normal::rv(&rng));
            
            std::map<Node*, std::vector<double>> simdat; //stores mean trait vectors
            std::vector<Node*> dpseq = t.getDownPassSequence();
            for (auto i = dpseq.rbegin(); i != dpseq.rend(); i++){
                Node* n = *i;
                if(n == t.getRoot()){
                    simdat.insert({n,rootMean});
                }else{
                    std::vector<double> ancMean =(simdat[n->getAncestor()]);
                    Eigen::MatrixXd var = sampledVCV * t.getBranchLength(n, n->getAncestor());
                    std::vector<double> draw = Probability::MultivariateNormal::rv(&rng, ancMean, &var);
                    simdat.insert({n,draw});
                }
            }
            
            std::vector<std::vector<double>> taxDat;
            std::vector<std::string> taxNames;
            for(auto a : simdat){
                if(a.first->getIsTip() == true){
                    taxNames.push_back(a.first->getName());
                    taxDat.push_back(a.second);
                }
            }
            
            Eigen::MatrixXd dat = Utility::EigenUtils::vectorMatrix2Eigen(taxDat);

            Tree internalTree = t;
            MultivariateBrownianMotionV2 mvBM(taxNames, &dat);
            mvBM.setTree(&internalTree);
            
//            TicToc tic("1000000 lnls");
//            for(int i = 0; i < 1000000; i++)
//                mvBM.lnLikelihood();
//            tic.end();
            
            Mcmc mcmc(1000000, 1000, 1000, &mvBM);
            TicToc loop("250000 generations");
            mcmc.run();
            loop.end();
            //evaluate VCV:
            UserSettings& settings = UserSettings::userSettings();
            ReadTSV r (settings.getOutputFile() + std::to_string(omp_get_thread_num()) + ".tsv", false, true);
            Eigen::MatrixXd rDat = r.getEigenMat();
            std::vector<std::string> cn = r.getColnames();
            for(int x = 0; x < cn.size(); x++){
                std::string s = cn[x];
                if(s.substr(0, 3) == "vcv"){
                    size_t rowPos = s.find_first_of("0123456789");
                    size_t commaPos = s.find(',');
                    int i = std::stoi(s.substr(rowPos, commaPos - rowPos));
                    int j = std::stoi(s.substr(commaPos + 1));
                    auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.1);
    //                auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.1);
                    if(sampledVCV(i,j) < impInterval.second && sampledVCV(i,j) > impInterval.first){
                        vcvInCredInt[i][j]++;
                        in++;
                    }else{
                    }
                }
            }
            
            std::cout<<std::endl;
            
            std::string old_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + ".tsv";
            std::string new_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + "mvBM" + std::to_string(cyc) + ".tsv";
            std::rename(old_name.c_str(), new_name.c_str());
            old_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + ".tre";
            new_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + "mvBM" + std::to_string(cyc) + ".tre";
            std::rename(old_name.c_str(), new_name.c_str());
        }
        std::cout << std::endl;
        for(int i = 0; i < ntraits; i++)
            for(int j = 0; j < ntraits; j++)
                std::cout << "(" + std::to_string(i) + "," + std::to_string(j) << ") " << vcvInCredInt[i][j] << "/" << std::to_string(nreps) << " | ";
        std::cout << std::endl;
        int total = 0;
        for(int i = 0; i < ntraits; i++)
            for(int j = 0; j < ntraits; j++)
                total += vcvInCredInt[i][j];
        std::cout << "total: " << total << " / " << nreps * ntraits * ntraits << " | " << (double)total / (double)(nreps * ntraits * ntraits) << std::endl;
    #endif
    
    #if 0
    // Multivariate Brownian Motion on LDDMM
    int nreps = 20;
    #pragma omp parallel for num_threads(10)
        for(int cyc = 0; cyc < nreps; cyc++){
            Tree t = Tree("(((Gorilla_beringei:2.558516,Gorilla_gorilla:2.558516):6.093717,((Homo_sapiens:0.568721,Neanderthal:0.538721):5.607159,(Pan_paniscus:2.333553,Pan_troglodytes:2.333553):3.842326):2.476353):6.480222,(Pongo_abelii:3.825854,Pongo_pygmaeus:3.825854):11.306601);");
            t.forceBinary();
                
            Eigen::MatrixXd rootShape2D = Utility::Shapes::generateUnitCirclePoints(5);
            SimShapeLDDMM s = SimShapeLDDMM(rootShape2D, &t);
            s.runSimulation(0.05, 0.1);
            
            std::map<Node*, Eigen::MatrixXd*> tipShapes = s.getTipShapes();
            std::vector<std::string> taxNames = s.getTipNames();
            Eigen::MatrixXd tipData = Eigen::MatrixXd::Zero(tipShapes.size(), 5 * 2);
            int idx = 0;
            for(auto a : tipShapes){
                Eigen::Map<const Eigen::RowVectorXd> flatRow(a.second->data(), a.second->size());
                tipData.row(idx) = flatRow;
                idx++;
            }
            MultivariateBrownianMotion mvBM(taxNames, &tipData);
            mvBM.lnLikelihood();
            Mcmc mcmc(1000000, 100, 100, &mvBM);
            mcmc.run();
            std::string old_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + ".tsv";
            std::string new_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + "LDDMMmvBM" + std::to_string(cyc) + ".tsv";
            std::rename(old_name.c_str(), new_name.c_str());
            old_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + ".tre";
            new_name = "/Users/levir/Documents/GitHub/phylo-master/output/out" + std::to_string(omp_get_thread_num()) + "LDDMMmvBM" + std::to_string(cyc) + ".tre";
            std::rename(old_name.c_str(), new_name.c_str());
        }
    #endif

    
    //PerikymataHSP v4
    # if 1
//    ReadCSV r = ReadCSV("/Users/levir/Documents/GitHub/Raskin_et_al_perikymata_hsp/LCdec3_10.csv", true, true);
//    ReadCSV r = ReadCSV("/Users/levir/Documents/GitHub/Raskin_et_al_perikymata_hsp/UI2dec3_10.csv", true, true);
    ReadCSV r = ReadCSV("/Users/levir/Documents/GitHub/Raskin_et_al_perikymata_hsp/UI2dec3_10_no_pongo.csv", true, true);
    std::vector<std::string> rawReadDatNames = r.getRownames();
    //remove quotation marks
    for(int i = 0 ; i < rawReadDatNames.size(); i++){
        std::string trimmed_name = rawReadDatNames[i];
        trimmed_name.erase(0, trimmed_name.find_first_not_of('"'));
        trimmed_name.erase(trimmed_name.find_last_not_of('"') + 1);
        rawReadDatNames[i] =trimmed_name;
    }
    Eigen::MatrixXd readDat = r.getEigenMat();
    Tree gatree = Tree("(((Gorilla_beringei:2.558516,Gorilla_gorilla:2.558516):6.093717,((Homo_sapiens:0.568721,Neanderthal:0.538721):5.607159,(Pan_paniscus:2.333553,Pan_troglodytes:2.333553):3.842326):2.476353):6.480222,(Pongo_abelii:3.825854,Pongo_pygmaeus:3.825854):11.306601);");
    int numChains = 10;
    std::vector<PhylogeneticModel*> perikymataModels;
    perikymataModels.resize(numChains);
    for(int i = 0; i < 10; i++)
        perikymataModels[i] = new PerikymataHSPv4(&gatree, rawReadDatNames, &readDat);
//        hsp.setVarianceCovarianceMatrix(sampledVCV);
    unsigned long ng = 100000000;
    MetropolisCoupledMcmc mcmc(ng, 1000, 1000, perikymataModels);
//    Mcmc mcmc(1000000, 1000, 1000, perikymataModels[0]);
    TicToc loop("1,000,000 generations");
    mcmc.run();
    loop.end();
    #endif
    
    #if 0
    //coverage check
    int ntraits = 10;
    int ntips = 8;
    int nIndividuals = 10;
    int nreps = 20;
    
    // Coverage tracking arrays
    int vcvInCredInt[ntraits][ntraits];
    int tipMeanInCredInt[ntips][ntraits];
    int tipVCVInCredInt[ntips][ntraits][ntraits];
    int imputedInCredInt = 0;

    // Initialize all coverage counters to zero
    for(int i = 0; i < ntraits; i++) {
        for(int j = 0; j < ntraits; j++) {
            vcvInCredInt[i][j] = 0;
            for(int k = 0; k < ntips; k++) {
                tipVCVInCredInt[k][i][j] = 0;
            }
        }
    }
    for(int k = 0; k < ntips; k++) {
        for(int i = 0; i < ntraits; i++) {
            tipMeanInCredInt[k][i] = 0;
        }
    }

    #pragma omp parallel for num_threads(10) shared(vcvInCredInt, tipMeanInCredInt, tipVCVInCredInt, imputedInCredInt)
    for(int cyc = 0; cyc < nreps; cyc++){
        //==============DATA==============//
        Tree rtree = Tree(ntips, 10.0);
            
        RandomVariable& rng = RandomVariable::randomVariableInstance();

        //mean evolution:
//        Eigen::VectorXd relRate = Eigen::VectorXd::Zero(ntraits);
//        Eigen::VectorXd dirichletA = Eigen::VectorXd::Constant(ntraits, 1.0);
//        Probability::Dirichlet::rv(&rng, dirichletA, relRate);
//        std::mt19937 mersenne;
//        LKJOnionSampler sampler(&mersenne, (size_t)ntraits);
//        Eigen::MatrixXd correlation = sampler.sampleEigen(1.0);
//        double averageRate = Probability::Normal::rv(&rng);
//        Eigen::VectorXd tmp = exp(averageRate) * relRate;
//        Eigen::MatrixXd sampledVCV = tmp.asDiagonal() * correlation * tmp.asDiagonal();
        Eigen::MatrixXd psi = Eigen::MatrixXd::Identity(ntraits, ntraits);;
        Eigen::MatrixXd sampledVCV = psi;
        Probability::InverseWishart::rv(&rng, sampledVCV, psi, ntraits+2);

        #pragma omp critical
        {
            std::cout << "========================================" << "\n";
            std::cout << "Evolutionary variance covariance matrix for thread " << omp_get_thread_num() << "\n";
            Utility::EigenUtils::printEigen(sampledVCV);
            std::cout << "========================================" << "\n";
        }
        
        std::vector<double> rootMean;
        for(int i = 0; i < sampledVCV.rows(); i++)
            rootMean.push_back(Probability::Normal::rv(&rng));
        std::map<Node*, std::vector<double>> simdat; //stores mean trait vectors
        std::vector<Node*> dpseq = rtree.getDownPassSequence();
        for (auto i = dpseq.rbegin(); i != dpseq.rend(); i++){
            Node* n = *i;
            if(n == rtree.getRoot()){
                simdat.insert({n,rootMean});
            }else{
                std::vector<double> ancMean =(simdat[n->getAncestor()]);
                Eigen::MatrixXd var = sampledVCV * rtree.getBranchLength(n, n->getAncestor());
                std::vector<double> draw = Probability::MultivariateNormal::rv(&rng, ancMean, &var);
                simdat.insert({n,draw});
            }
        }

        std::vector<std::vector<double>> taxDat;
        std::vector<std::string> taxNames;
        for(auto a : simdat){
            if(a.first->getIsTip() == true){
                taxNames.push_back(a.first->getName());
                taxDat.push_back(a.second);
            }
        }
        
        //vector of tip means
        Eigen::MatrixXd dat = Utility::EigenUtils::vectorMatrix2Eigen(taxDat);
        
        // Store true intraspecific parameters
        std::map<std::string, Eigen::VectorXd> trueTipMeans;
        std::map<std::string, Eigen::MatrixXd> trueTipVCVs;
        std::map<std::string, int> tipNameToIndex;
        
        //sample individuals at each tip
        Eigen::MatrixXd indTipDat = Eigen::MatrixXd::Zero(ntips*nIndividuals, ntraits);
        int idx = 0;
        std::vector<std::string> rowTipNames;
        for(int i = 0; i < ntips; i++){
            //sample tip VCV
            double dof = ntraits + 2; //such that the mean is the scale matrix
            double lambda = dof - ntraits - 1;
            Eigen::VectorXd mu0 = Eigen::VectorXd::Zero(ntraits);
            Eigen::MatrixXd psi = Eigen::MatrixXd::Constant(ntraits, ntraits, 1e-6);
            psi.diagonal().array() = 1.0;
            psi = Eigen::VectorXd::Ones(ntraits).asDiagonal() * psi * Eigen::VectorXd::Ones(ntraits).asDiagonal();
            Eigen::MatrixXd tipVCV = Probability::InverseWishart::rv(&rng, psi, dof);
            
            std::string tipName = taxNames[i];
            Eigen::VectorXd tipMean = dat.row(i);
            
            // Store true values
            trueTipMeans[tipName] = tipMean;
            trueTipVCVs[tipName] = tipVCV;
            tipNameToIndex[tipName] = i;
            
            for(int j = 0; j < nIndividuals; j++){
                Eigen::VectorXd tmprow = Probability::MultivariateNormal::rv(&rng, tipMean, &tipVCV);
                indTipDat.row(idx) = tmprow;
                rowTipNames.push_back(tipName);
                idx++;
            }
        }
        
        //add in one NaN
        int missingI = (int)(rng.uniformRv() * indTipDat.rows());
        int missingJ = (int)(rng.uniformRv() * indTipDat.cols());
        double trueVal = indTipDat(missingI, missingJ);
        indTipDat(missingI, missingJ) = std::numeric_limits<double>::quiet_NaN();
        
        //==============MCMC==============//
        std::vector<PhylogeneticModel*> perikymataModels;
        perikymataModels.resize(10);
        for(int i = 0; i < 10; i++)
            perikymataModels[i] = new PerikymataHSPv4(&rtree, rowTipNames, &indTipDat);
        MetropolisCoupledMcmc mcmc(10000000, 1000, 1000, perikymataModels);
        TicToc t("optims");
//        Mcmc mcmc(10000000, 1000, 1000, perikymataModels[0]);
        mcmc.run();
        t.end();
        
        //==============COVERAGE EVALUATION==============//
        UserSettings& settings = UserSettings::userSettings();
        ReadTSV r(settings.getOutputFile() + std::to_string(omp_get_thread_num()) + ".tsv", false, true);
        Eigen::MatrixXd rDat = r.getEigenMat();
        std::vector<std::string> cn = r.getColnames();
        
        // Evaluate evolutionary VCV coverage
        int evoVCVCovered = 0;
        for(int x = 0; x < cn.size(); x++){
            std::string s = cn[x];
            if(s.substr(0, 3) == "vcv"){
                size_t rowPos = s.find_first_of("0123456789");
                size_t commaPos = s.find(',');
                int i = std::stoi(s.substr(rowPos, commaPos - rowPos));
                int j = std::stoi(s.substr(commaPos + 1));
                auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.1);
                if(sampledVCV(i,j) < impInterval.second && sampledVCV(i,j) > impInterval.first){
                    #pragma omp atomic
                    vcvInCredInt[i][j]++;
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
                std::string paramName = "intraspecificMean" + tipName + std::to_string(traitIdx);
                
                auto it = std::find(cn.begin(), cn.end(), paramName);
                if(it != cn.end()) {
                    int colIdx = std::distance(cn.begin(), it);
                    auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(colIdx), 0.1);
                    
                    if(trueMean(traitIdx) < impInterval.second && trueMean(traitIdx) > impInterval.first) {
                        #pragma omp atomic
                        tipMeanInCredInt[tipIdx][traitIdx]++;
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
                    std::string paramName = "intraspecificVCV" + tipName + std::to_string(i) + "," + std::to_string(j);
                    
                    auto it = std::find(cn.begin(), cn.end(), paramName);
                    if(it != cn.end()) {
                        int colIdx = std::distance(cn.begin(), it);
                        auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(colIdx), 0.1);
                        
                        if(trueVCV(i,j) < impInterval.second && trueVCV(i,j) > impInterval.first) {
                            #pragma omp atomic
                            tipVCVInCredInt[tipIdx][i][j]++;
                            tipVCVCovered++;
                        }
                    }
                }
            }
        }
        
        // Evaluate missing data imputation
        for(int x = 0; x < cn.size(); x++){
            std::string s = cn[x];
            if(s.substr(0, 4) == "pkmm"){
                auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.1);
                if(trueVal < impInterval.second && trueVal > impInterval.first){
                    #pragma omp atomic
                    imputedInCredInt++;
                }
            }
        }
        
        // Rename output files
        std::string old_name = settings.getOutputFile() + std::to_string(omp_get_thread_num()) + ".tsv";
        std::string new_name = settings.getOutputFile() + std::to_string(omp_get_thread_num()) + "mvBM" + std::to_string(cyc) + ".tsv";
        std::rename(old_name.c_str(), new_name.c_str());
        old_name = settings.getOutputFile() + std::to_string(omp_get_thread_num()) + ".tre";
        new_name = settings.getOutputFile() + std::to_string(omp_get_thread_num()) + "mvBM" + std::to_string(cyc) + ".tre";
        std::rename(old_name.c_str(), new_name.c_str());
    }

    //==============FINAL SUMMARY==============//
    std::cout << "\n========================================" << std::endl;
    std::cout << "FINAL COVERAGE RESULTS" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Evolutionary VCV summary
    std::cout << "--- EVOLUTIONARY VCV MATRIX ---" << std::endl;
    for(int i = 0; i < ntraits; i++) {
        for(int j = 0; j < ntraits; j++) {
            std::cout << "vcv[" << i << "," << j << "]: "
                      << vcvInCredInt[i][j] << "/" << nreps
                      << " = " << (double)vcvInCredInt[i][j] / nreps << std::endl;
        }
    }
    int total = 0;
    for(int i = 0; i < ntraits; i++)
        for(int j = 0; j < ntraits; j++)
            total += vcvInCredInt[i][j];
    std::cout << "Overall Evo VCV Coverage: " << total << " / " << (nreps * ntraits * ntraits)
              << " = " << (double)total / (nreps * ntraits * ntraits) << "\n" << std::endl;

    // Tip means summary
    std::cout << "--- INTRASPECIFIC MEANS ---" << std::endl;
    int totalTipMeans = 0;
    for(int k = 0; k < ntips; k++) {
        int tipTotal = 0;
        for(int i = 0; i < ntraits; i++) {
            tipTotal += tipMeanInCredInt[k][i];
            totalTipMeans += tipMeanInCredInt[k][i];
        }
        std::cout << "Tip " << k << " mean coverage: " << tipTotal << "/" << (nreps * ntraits)
                  << " = " << (double)tipTotal / (nreps * ntraits) << std::endl;
    }
    std::cout << "Overall Tip Mean Coverage: " << totalTipMeans << " / " << (nreps * ntips * ntraits)
              << " = " << (double)totalTipMeans / (nreps * ntips * ntraits) << "\n" << std::endl;

    // Tip VCVs summary
    std::cout << "--- INTRASPECIFIC VCV MATRICES ---" << std::endl;
    int totalTipVCVs = 0;
    for(int k = 0; k < ntips; k++) {
        int tipTotal = 0;
        for(int i = 0; i < ntraits; i++) {
            for(int j = 0; j < ntraits; j++) {
                tipTotal += tipVCVInCredInt[k][i][j];
                totalTipVCVs += tipVCVInCredInt[k][i][j];
            }
        }
        std::cout << "Tip " << k << " VCV coverage: " << tipTotal << "/" << (nreps * ntraits * ntraits)
                  << " = " << (double)tipTotal / (nreps * ntraits * ntraits) << std::endl;
    }
    std::cout << "Overall Tip VCV Coverage: " << totalTipVCVs << " / " << (nreps * ntips * ntraits * ntraits)
              << " = " << (double)totalTipVCVs / (nreps * ntips * ntraits * ntraits) << "\n" << std::endl;

    // Missing data summary
    std::cout << "--- MISSING DATA IMPUTATION ---" << std::endl;
    std::cout << "Missing Data Coverage: " << imputedInCredInt << "/" << nreps
              << " = " << (double)imputedInCredInt / nreps << "\n" << std::endl;

    // Overall summary
    std::cout << "========================================" << std::endl;
    std::cout << "OVERALL SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;
    int totalParams = (nreps * ntraits * ntraits) +                    // Evo VCV
                      (nreps * ntips * ntraits) +                       // Tip means
                      (nreps * ntips * ntraits * ntraits) +            // Tip VCVs
                      nreps;                                            // Missing data
    int totalCovered = total + totalTipMeans + totalTipVCVs + imputedInCredInt;
    std::cout << "Total Parameters Tested: " << totalParams << std::endl;
    std::cout << "Total Parameters Covered: " << totalCovered << std::endl;
    std::cout << "Overall Coverage: " << (double)totalCovered / totalParams << std::endl;
    std::cout << "========================================\n" << std::endl;

    #endif
    
    //Hominin disparity
    #if 0
    
    // Simulation parameters:
    //constant parameters
    int numReps = 100;
    double lddmmSigma = 1.0;
    double betaDistAlpha = 1.0;
    double betaDistBeta = 1.0;
    int threads = 10;
    
    //changing parameters
    std::vector<double> alphas = {0.1, 0.2, 0.3, 0.4};
    std::vector<int> numAdditionalTaxa = {3, 5, 8, 11, 15, 20, 25, 50, 100, 144}; //abstract explicitly mentions 5, 11, 25,50, and 144
//    std::vector<int> numLandmarks = {10, 25, 50};
    std::vector<int> numLandmarks = {25, 50};
    
    std::string base ="/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/";
    std::string treeIn = "/Users/levir/Documents/GitHub/HomininTaxicDiversity/data/sampledTrees.tsv";

    int totalIter = numLandmarks.size() * alphas.size() * numAdditionalTaxa.size();
    int iter = 0;
    auto loopStart = std::chrono::steady_clock::now();

    for(int lm : numLandmarks){
        for(double a : alphas){
            for(int nt : numAdditionalTaxa){
                std::cout << "Current simulation:  " << "\n";
                std::cout << "Num landmarks:       " << lm << "\n";
                std::cout << "LDDMM alpha:         " << a << "\n";
                std::cout << "Num additional taxa: " << nt << "\n";
                
                // ETA code
                if(iter > 0){
                    auto now = std::chrono::steady_clock::now();
                    double elapsed = std::chrono::duration<double>(now - loopStart).count();
                    double avgPerIter = elapsed / iter;
                    double eta = avgPerIter * (totalIter - iter);
                    int etaMin = (int)(eta / 60);
                    int etaSec = (int)(eta) % 60;
                    std::cout << "[" << iter << "/" << totalIter << "] ETA: " << etaMin << "m " << etaSec << "s" << " | Average per iteration: " << avgPerIter  <<  "\n";
                }
                
                int numAddTaxa = nt;
                double lddmmAlpha = a;
                int numLM = lm;
                
                
                //helperlambda for tree out/shape out
                auto fmt_path = [](int lm, int nt, double a) {
                    char buf[256];
                    std::snprintf(buf, sizeof(buf), "numLandmarks%d/numAdditionalTaxa%d/alpha%.2f/", lm, nt, a);
                    return std::string(buf);
                };
                
                std::string treeOut =   base +
                                        "simulatedTrees/" +
                                        fmt_path(lm , nt , a);
                std::string shapeOut =  base +
                                        "simulatedShapes/"+
                                        fmt_path(lm , nt , a);
                                        
                std::error_code ec;
                std::filesystem::create_directories(treeOut, ec);
                if(ec) std::cerr << "Failed to create treeOut: " << ec.message() << "\n";
                std::filesystem::create_directories(shapeOut, ec);
                if(ec) std::cerr << "Failed to create shapeOut: " << ec.message() << "\n";
                
                HomininDisparityAnalysis disp(treeIn, treeOut, shapeOut);
                disp.run(numAddTaxa, numReps, lddmmSigma, lddmmAlpha, betaDistAlpha, betaDistBeta, numLM, threads);
                iter++;
            }
        }
    }
    
    #endif
    
    return 0;
}
