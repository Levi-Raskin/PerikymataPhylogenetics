#include "Mcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "TicToc.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "WriteTSV.hpp"

#include <iomanip>
#include <iostream>
#include </usr/local/include/omp.h>


Mcmc::Mcmc(int ng, int pf, int sf, PhylogeneticModel* m) : numCycles(ng), printFrequency(pf), sampleFrequency(sf), model(m), treeStrm(nullptr) {
    
}

Mcmc::~Mcmc(void) {

    if (treeStrm != nullptr)
        delete treeStrm;
}

void Mcmc::run(void) {
    //phylogenetic model
    if(model != nullptr){
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        
        openFiles();

        double curLnL = model->lnLikelihood();
        double curLnP = model->lnPriorProbability();

        if(numCycles >= std::numeric_limits<unsigned long>::max())
            Msg::error("numCycles requested in greater than largest possible value");

        for (unsigned long n=1; n<=numCycles; n++)
            {
//            TicToc p1 ("update");
            double lnProposalRatio = model->update();
//            p1.end();
//            TicToc p2 ("lnl");
            double newLnL = model->lnLikelihood();
//            p2.end();
            
            //calc power here
            double lnLikelihoodRatio = newLnL - curLnL;
//            TicToc p3 ("lnp");
            double newLnP = model->lnPriorProbability();
//            p3.end();
            double lnPriorRatio = newLnP - curLnP;
            double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

//            std::cout << "lnLikelihoodRatio: " << lnLikelihoodRatio << std::endl;
//            std::cout << "lnPriorRatio: " << lnPriorRatio << std::endl;
//            std::cout << "lnProposalRatio: " << lnProposalRatio << std::endl;
//            std::cout << "lnR: " << lnR << std::endl;
            
            // accept or reject the proposed state
            bool acceptMove = false;
            if (log(rng.uniformRv()) < lnR)
                acceptMove = true;
                
            if (n % printFrequency == 0)
                {
                std::cout << std::fixed << std::setprecision(2);
                std::cout << n << " -- " << curLnL << " -> " << newLnL;
                model->print();
                }
                
            // adjust the chain accordingly
            if (acceptMove == true)
                {
                curLnL = newLnL;
                curLnP = newLnP;
                model->updateForAcceptance();
                }
            else
                {
                model->updateForRejection();
                }
                
                //choose two chains and swap
                
            if (n == 1 || n == numCycles || n % sampleFrequency == 0 )
                sample(n, curLnL);
            }
    }else{
        Msg::error("model somehow not instantiated");
    }
}

void Mcmc::openFiles(void) {

    UserSettings& settings = UserSettings::userSettings();
    std::string treeFileName = settings.getOutputFile();
//    std::string parmFileName = settings.getOutputFile();
    treeFileName += std::to_string(omp_get_thread_num());
    treeFileName += ".tre";
//    parmFileName += ".csv";
    
    treeStrm = new std::ofstream(treeFileName.c_str(), std::ios::out);
//    parmStrm = new std::ofstream(parmFileName.c_str(), std::ios::out);
    if (treeStrm == nullptr)
        Msg::error("Could not open tree and/or parameter files");
}

void Mcmc::sample(unsigned long n, double lnL, double lnP) {
    UserSettings& settings = UserSettings::userSettings();
    std::string tracerFileName = settings.getOutputFile();
    std::string parmFileName = tracerFileName;
//    treeFileName += ".tre";
    parmFileName += std::to_string(omp_get_thread_num());
    parmFileName += ".tsv";
    WriteTSV w = WriteTSV();
    
    if(n == 1){
        w.addFilepath(parmFileName, true);
        std::vector<std::string> cn = {"n", "lnL", "lnP"};
        std::vector<std::string> headStr = model->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        w.addColumnNamesTSV(cn);
        *treeStrm << "#NEXUS\n\n";
        *treeStrm << "begin trees;\n";
    }
    else
        w.addFilepath(parmFileName, false);
        
    std::vector<double> dat = {(double)n, lnL, lnP};
    std::vector<double> parmStr = model->getParameterString();
    dat.insert( dat.end(), parmStr.begin(), parmStr.end() );
    w.appendDataTSV(dat);

    Tree* t = model->getTree();
    std::string treeStr = t->getNewickString();
    *treeStrm << "   tree tree_" << n << " = " << treeStr << "\n";

    //    std::vector<double> dat = {(double)n, lnL, parmValues[0]};
    //    w.appendDataTSV(dat);

    if (n == numCycles)
        {
        *treeStrm << "end;\n";
        w.closeTSV();
        }
    /*
    UserSettings& settings = UserSettings::userSettings();
    std::string tracerFileName = settings.getOutputFile();
//    tracerFileName += ".tsv";
//    WriteTSV w = WriteTSV();
//    if(n == 1){
//        w.addFilepath(tracerFileName, true);
//        w.addColumnNamesTSV({"iteration", "lh", "k0"});
//    }
//    else
//        w.addFilepath(tracerFileName, false);

    if (n == 1)
        {
        *treeStrm << "#NEXUS\n\n";
        *treeStrm << "begin trees;\n";
        
        std::vector<std::string> header = model->getHeaderString();
        *parmStrm << "N\tlnL\t";
        for (int i=0; i<header.size(); i++)
            *parmStrm << header[i] << '\t';
        *parmStrm << '\n';
        }
        
    Tree* t = model->getTree();
    std::string treeStr = t->getNewickString();
    *treeStrm << "   tree tree_" << n << " = " << treeStr << "\n";
    
    std::vector<double> parmValues = model->getParameterString();
    *parmStrm << n << '\t' << lnL << '\t';
    for (int i=0; i<parmValues.size(); i++)
        *parmStrm << parmValues[i] << '\t';
    *parmStrm << '\n';
    
//    std::vector<double> dat = {(double)n, lnL, parmValues[0]};
//    w.appendDataTSV(dat);
        
    if (n == numCycles)
        {
        *treeStrm << "end;\n";
//        w.closeTSV();
        }
    */
}

void Mcmc::sample(unsigned long n, double lnL) {
    UserSettings& settings = UserSettings::userSettings();
    std::string tracerFileName = settings.getOutputFile();
    std::string parmFileName = tracerFileName;
//    treeFileName += ".tre";
    parmFileName += std::to_string(omp_get_thread_num());
    parmFileName += ".tsv";
    WriteTSV w = WriteTSV();
    
    if(n == 1){
        w.addFilepath(parmFileName, true);
        std::vector<std::string> cn = {"n", "lnL"};
        std::vector<std::string> headStr = model->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        w.addColumnNamesTSV(cn);
        *treeStrm << "#NEXUS\n\n";
        *treeStrm << "begin trees;\n";
    }
    else
        w.addFilepath(parmFileName, false);
        
    std::vector<double> dat = {(double)n, lnL};
    std::vector<double> parmStr = model->getParameterString();
    dat.insert( dat.end(), parmStr.begin(), parmStr.end() );
    w.appendDataTSV(dat);

    Tree* t = model->getTree();
    if(t != nullptr){
        std::string treeStr = t->getNewickString();
        *treeStrm << "   tree tree_" << n << " = " << treeStr << "\n";
    }
    //    std::vector<double> dat = {(double)n, lnL, parmValues[0]};
    //    w.appendDataTSV(dat);

    if (n == numCycles)
        {
        *treeStrm << "end;\n";
        w.closeTSV();
        }
    /*
    UserSettings& settings = UserSettings::userSettings();
    std::string tracerFileName = settings.getOutputFile();
//    tracerFileName += ".tsv";
//    WriteTSV w = WriteTSV();
//    if(n == 1){
//        w.addFilepath(tracerFileName, true);
//        w.addColumnNamesTSV({"iteration", "lh", "k0"});
//    }
//    else
//        w.addFilepath(tracerFileName, false);

    if (n == 1)
        {
        *treeStrm << "#NEXUS\n\n";
        *treeStrm << "begin trees;\n";
        
        std::vector<std::string> header = model->getHeaderString();
        *parmStrm << "N\tlnL\t";
        for (int i=0; i<header.size(); i++)
            *parmStrm << header[i] << '\t';
        *parmStrm << '\n';
        }
        
    Tree* t = model->getTree();
    std::string treeStr = t->getNewickString();
    *treeStrm << "   tree tree_" << n << " = " << treeStr << "\n";
    
    std::vector<double> parmValues = model->getParameterString();
    *parmStrm << n << '\t' << lnL << '\t';
    for (int i=0; i<parmValues.size(); i++)
        *parmStrm << parmValues[i] << '\t';
    *parmStrm << '\n';
    
//    std::vector<double> dat = {(double)n, lnL, parmValues[0]};
//    w.appendDataTSV(dat);
        
    if (n == numCycles)
        {
        *treeStrm << "end;\n";
//        w.closeTSV();
        }
    */
}
