#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "ReadTSV.hpp"
#include "SimulateData.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "Utility.hpp"

#include <cmath>
#include <iostream>
#include <fstream>

SimulateData::SimulateData(void) :
    tree(nullptr),
    incrementedElapsed(0.0),
    incrementRunning(false),
    rankFileInitialized(false),
    rankBurninFraction(0.1),
    rankThinStride(1),
    trials(0){
    
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

    evoVCVRepCoverage.reserve(nreps);
    tipVCVRepCoverage.reserve(nreps);
    tipMeanRepCoverage.reserve(nreps);
    missingRepCoverage.reserve(nreps);
    totalRepCoverage.reserve(nreps);
    
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
    Probability::InverseWishart::rv(&rng, sampledEvoVCV, psi.llt().matrixL(), priorDOF);
    
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

void SimulateData::startIncrement(void){
    incrementStart = std::chrono::steady_clock::now();
    incrementRunning = true;
}

void SimulateData::endIncrement(void){
    if(incrementRunning == false){
        Msg::warning("endIncrement() called without a matching startIncrement(); ignoring");
        return;
    }
    auto now = std::chrono::steady_clock::now();
    incrementedElapsed += std::chrono::duration<double>(now - incrementStart).count();
    incrementRunning = false;
}

Eigen::VectorXd SimulateData::postBurninThinned(const Eigen::VectorXd& col) const {
    int n = (int)col.size();
    int discard = (int)(n * rankBurninFraction);
    if(discard >= n)
        Msg::error("postBurninThinned: rankBurninFraction leaves no post-burn-in samples");

    std::vector<double> kept;
    kept.reserve((n - discard) / rankThinStride + 1);
    for(int i = discard; i < n; i += rankThinStride)
        kept.push_back(col(i));

    Eigen::VectorXd out(kept.size());
    for(size_t i = 0; i < kept.size(); i++)
        out(i) = kept[i];
    return out;
}

Eigen::VectorXd SimulateData::vcvTraceTrace(const Eigen::MatrixXd& rDat, const std::vector<std::string>& cn,
                                              const std::string& vcvPrefix, int nTraits) const {
    Eigen::VectorXd result;
    for(int i = 0; i < nTraits; i++){
        std::string colName = vcvPrefix + "_(" + std::to_string(i) + "," + std::to_string(i) + ")";
        auto it = std::find(cn.begin(), cn.end(), colName);
        if(it == cn.end())
            Msg::error("vcvTraceTrace: could not find diagonal column " + colName);
        int colIdx = (int)std::distance(cn.begin(), it);
        Eigen::VectorXd trimmed = postBurninThinned(rDat.col(colIdx));
        if(i == 0) result = trimmed;
        else       result += trimmed;
    }
    return result;
}

int SimulateData::computeRank(const Eigen::VectorXd& posterior, double trueVal) const {
    int rank = 0;
    for(int l = 0; l < posterior.size(); l++)
        if(posterior(l) < trueVal)
            rank++;
    return rank;
}

void SimulateData::writeRankRow(const std::string& label, int rank, int L){
    double normalizedRank = (double)rank / (double)(L + 1);
    rankOut << trials << "," << label << "," << rank << "," << L << "," << normalizedRank << "\n";
}

void SimulateData::checkRankUniformity(const Eigen::MatrixXd& rDat, const std::vector<std::string>& cn){
    UserSettings& settings = UserSettings::userSettings();

    if(rankFileInitialized == false){
        rankOut.open(settings.getOutputFile() + "rank_uniform_test.csv");
        if(!rankOut.is_open())
            Msg::error("Could not open rank uniformity output file");
        rankOut << "replicate_id,parameter_label,rank,L,normalized_rank\n";
        rankFileInitialized = true;
    }

    // --- evolutionary VCV: rank of the TRUE trace among posterior trace draws ---
    if(settings.getWithPhylogeny() == true){
        Eigen::VectorXd evoTrace = vcvTraceTrace(rDat, cn, "evo_vcv", ntraits);
        double trueTrace = sampledEvoVCV.trace();
        int rank = computeRank(evoTrace, trueTrace);
        writeRankRow("evo_vcv_trace", rank, (int)evoTrace.size());
    }

    // --- per-taxon intraspecific VCV: rank of TRUE trace among posterior trace draws ---
    for(auto& tipEntry : trueTipVCVs){
        const std::string& tipName = tipEntry.first;
        const Eigen::MatrixXd& trueVCV = tipEntry.second;

        std::string firstCol = tipName + "_vcv_(0,0)";
        if(std::find(cn.begin(), cn.end(), firstCol) == cn.end())
            continue;   // e.g. a singleton-observation tip has no VCV parameter; skip cleanly

        Eigen::VectorXd tipTrace = vcvTraceTrace(rDat, cn, tipName + "_vcv", ntraits);
        double trueTrace = trueVCV.trace();
        int rank = computeRank(tipTrace, trueTrace);
        writeRankRow(tipName + "_vcv_trace", rank, (int)tipTrace.size());
    }

    // --- per-taxon, per-trait intraspecific mean: rank for EACH component ---
    for(auto& tipEntry : trueTipMeans){
        const std::string& tipName = tipEntry.first;
        const Eigen::VectorXd& trueMean = tipEntry.second;

        std::string firstCol = tipName + "_mean_0";
        if(std::find(cn.begin(), cn.end(), firstCol) == cn.end())
            continue;   // singleton-observation tip: mean is fixed, not sampled

        for(int k = 0; k < ntraits; k++){
            std::string colName = tipName + "_mean_" + std::to_string(k);
            auto it = std::find(cn.begin(), cn.end(), colName);
            if(it == cn.end()){
                Msg::warning("checkRankUniformity: could not find column " + colName + "; skipping");
                continue;
            }
            int colIdx = (int)std::distance(cn.begin(), it);
            Eigen::VectorXd posterior = postBurninThinned(rDat.col(colIdx));
            int rank = computeRank(posterior, trueMean(k));
            writeRankRow(tipName + "_mean_" + std::to_string(k), rank, (int)posterior.size());
        }
    }

    rankOut.flush();
}

void SimulateData::checkCredInt(void){
    UserSettings& settings = UserSettings::userSettings();
    ReadTSV r(settings.getOutputFile() + "Outfile.tsv", false, true);
    Eigen::MatrixXd rDat = r.getEigenMat();
    std::vector<std::string> cn = r.getColnames();
    
    checkRankUniformity(rDat, cn);
    
    // Evaluate evolutionary VCV coverage
    int evoVCVCovered = 0;
    for(int x = 0; x < cn.size(); x++){
        std::string s = cn[x];
        if(s.substr(0, 7) == "evo_vcv"){
            size_t rowPos = s.find_first_of("0123456789");
            size_t commaPos = s.find(',');
            int i = std::stoi(s.substr(rowPos, commaPos - rowPos));
            int j = std::stoi(s.substr(commaPos + 1));
            auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.5);
            if(sampledEvoVCV(i,j) < impInterval.second && sampledEvoVCV(i,j) > impInterval.first){
                vcvInCredInt(i,j)++;
                evoVCVCovered++;
            }
        }
    }
    // Evaluate tip mean coverage
    int tipMeanCoveredThisRep = 0;
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
                auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(colIdx), 0.5);
                
                if(trueMean(traitIdx) < impInterval.second && trueMean(traitIdx) > impInterval.first) {
                    tipMeanInCredInt(tipIdx,traitIdx)++;
                    tipMeanCovered++;
                }
            }
        }
        tipMeanCoveredThisRep += tipMeanCovered;
    }
    
    // Evaluate tip VCV coverage
    int tipVCVCoveredThisRep = 0;
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
                    auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(colIdx), 0.5);
                    
                    if(trueVCV(i,j) < impInterval.second && trueVCV(i,j) > impInterval.first) {
                        tipVCVInCredInt[tipIdx](i,j)++;
                        tipVCVCovered++;
                    }
                }
            }
        }
        tipVCVCoveredThisRep += tipVCVCovered;
    }
    
    // Evaluate missing data imputation
    int missingCoveredThisRep = 0;
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
                auto impInterval = Utility::Bayesian::credibleIntervalBurnIn(rDat.col(x), 0.5);
                if(trueVal > impInterval.first && trueVal < impInterval.second){
                    int impIdx = (int)std::distance(trueMissingValues.begin(), it);
                    imputedInCredInt(impIdx)++;
                    missingCoveredThisRep++;
                }
            }
        }
    }

    double evoRateThisRep     = (ntraits > 0) ? (double)evoVCVCovered / (double)(ntraits * ntraits) : std::nan("");
    double tipVCVRateThisRep  = (ntips > 0 && ntraits > 0) ? (double)tipVCVCoveredThisRep / (double)(ntips * ntraits * ntraits) : std::nan("");
    double tipMeanRateThisRep = (ntips > 0 && ntraits > 0) ? (double)tipMeanCoveredThisRep / (double)(ntips * ntraits) : std::nan("");

    evoVCVRepCoverage.push_back(evoRateThisRep);
    tipVCVRepCoverage.push_back(tipVCVRateThisRep);
    tipMeanRepCoverage.push_back(tipMeanRateThisRep);

    if(nimp > 0){
        double missingRateThisRep = (double)missingCoveredThisRep / (double)nimp;
        missingRepCoverage.push_back(missingRateThisRep);
    }

    int numThisRep = 0, denThisRep = 0;
    if(settings.getWithIntraspecific() == false){
        numThisRep = evoVCVCovered;
        denThisRep = ntraits * ntraits;
    } else if(settings.getWithPhylogeny() == false){
        numThisRep = tipVCVCoveredThisRep + tipMeanCoveredThisRep + missingCoveredThisRep;
        denThisRep = (ntips * ntraits * ntraits) + (ntips * ntraits) + nimp;
    } else {
        numThisRep = evoVCVCovered + tipVCVCoveredThisRep + tipMeanCoveredThisRep + missingCoveredThisRep;
        denThisRep = (ntraits * ntraits) + (ntips * ntraits * ntraits) + (ntips * ntraits) + nimp;
    }
    double totalRateThisRep = (denThisRep > 0) ? (double)numThisRep / (double)denThisRep : std::nan("");
    totalRepCoverage.push_back(totalRateThisRep);
}

std::pair<double,double> SimulateData::meanAndSE(const std::vector<double>& v) const {
    std::vector<double> clean;
    clean.reserve(v.size());
    for(double x : v)
        if(!std::isnan(x))
            clean.push_back(x);

    if(clean.empty())
        return std::make_pair(std::nan(""), std::nan(""));

    double sum = 0.0;
    for(double x : clean)
        sum += x;
    double mean = sum / (double)clean.size();

    if(clean.size() < 2)
        return std::make_pair(mean, std::nan(""));

    double ss = 0.0;
    for(double x : clean)
        ss += (x - mean) * (x - mean);
    double sd = std::sqrt(ss / (double)(clean.size() - 1));
    double se = sd / std::sqrt((double)clean.size());
    
    return std::make_pair(mean, se);
}

void SimulateData::print(void){
    int total = 0;
    for(auto& m : tipVCVInCredInt)
        total += m.sum();

    std::pair<double,double> evoStats     = meanAndSE(evoVCVRepCoverage);
    std::pair<double,double> tipVCVStats  = meanAndSE(tipVCVRepCoverage);
    std::pair<double,double> tipMeanStats = meanAndSE(tipMeanRepCoverage);
    std::pair<double,double> missingStats = meanAndSE(missingRepCoverage);
    std::pair<double,double> totalStats   = meanAndSE(totalRepCoverage);

    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Pooled coverage (all element x replicate checks)" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Evolutionary VCV coverage:            " << vcvInCredInt.sum() << "/" << (trials * ntraits * ntraits) << "\t\t | (" << (double)vcvInCredInt.sum() / (trials * ntraits * ntraits) << ")" << "\n";
    std::cout << "Tip VCV coverage:                     " << total << "/" << (trials * ntips * ntraits * ntraits) << "\t\t | (" << (double)total / (trials * ntips * ntraits * ntraits) << ")" << "\n";
    std::cout << "Tip mean coverage:                    " << tipMeanInCredInt.sum() << "/" << (trials * ntips * ntraits) << "\t\t | (" << (double)tipMeanInCredInt.sum() / (trials * ntips * ntraits) << ")" << "\n";
    std::cout << "Missing data coverage:                " << imputedInCredInt.sum() << "/" << (trials * nimp) << "\t\t | (" << (double)imputedInCredInt.sum() / (trials * nimp) << ")" << "\n";
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Per-replicate coverage (mean +/- SE across " << trials << " replicates, SE = SD/sqrt(n reps));" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Total coverage:                       " << totalStats.first << " +/- " << totalStats.second << "\n";
    std::cout << "Evolutionary VCV coverage:            " << evoStats.first << " +/- " << evoStats.second << "\n";
    std::cout << "Tip VCV coverage:                     " << tipVCVStats.first << " +/- " << tipVCVStats.second << "\n";
    std::cout << "Tip mean coverage:                    " << tipMeanStats.first << " +/- " << tipMeanStats.second << "\n";
    std::cout << "Missing data coverage:                " << missingStats.first << " +/- " << missingStats.second << "\n";
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

    std::pair<double,double> evoStats     = meanAndSE(evoVCVRepCoverage);
    std::pair<double,double> tipVCVStats  = meanAndSE(tipVCVRepCoverage);
    std::pair<double,double> tipMeanStats = meanAndSE(tipMeanRepCoverage);
    std::pair<double,double> missingStats = meanAndSE(missingRepCoverage);
    std::pair<double,double> totalStats   = meanAndSE(totalRepCoverage);

    log << "-----------------------------------------------------------------------\n";
    log << "Per-replicate coverage (mean +/- SE across " << trials << " replicates; SE = SD/sqrt(n reps))\n";
    log << "-----------------------------------------------------------------------\n";
    log << "Total coverage: " << totalStats.first << " +/- " << totalStats.second << "\n";
    log << "-----------------------------------------------------------------------\n";
    log << "Evolutionary VCV coverage:            " << evoStats.first << " +/- " << evoStats.second << "\n";
    log << "Tip VCV coverage:                     " << tipVCVStats.first << " +/- " << tipVCVStats.second << "\n";
    log << "Tip mean coverage:                    " << tipMeanStats.first << " +/- " << tipMeanStats.second << "\n";
    log << "Missing data coverage:                " << missingStats.first << " +/- " << missingStats.second << "\n";

    log.close();
    
    if(rankFileInitialized)
        rankOut.close();
}
