#include "Mcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "WriteTSV.hpp"

#include <iomanip>
#include <iostream>


Mcmc::Mcmc(int ng, int pf, int sf, PhylogeneticModel* m, std::string of) : numCycles(ng), printFrequency(pf), sampleFrequency(sf), model(m), treeStrm(nullptr), outfile(of) {
    
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

        for (int n=1; n<=numCycles; n++)
            {
            double lnProposalRatio = model->update();
            double newLnL = model->lnLikelihood();
            double lnLikelihoodRatio = newLnL - curLnL;
            double newLnP = model->lnPriorProbability();
            double lnPriorRatio = newLnP - curLnP;
            double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

            std::cout << "lnLikelihoodRatio: " << lnLikelihoodRatio << std::endl;
            std::cout << "lnPriorRatio: " << lnPriorRatio << std::endl;
            std::cout << "lnProposalRatio: " << lnProposalRatio << std::endl;
            std::cout << "lnR: " << lnR << std::endl;
            
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
                
            if (n == 1 || n == numCycles || n % sampleFrequency == 0 )
                sample(n, curLnL);
            }
    }else{
        Msg::error("model somehow not instantiated");
    }
}

void Mcmc::openFiles(void) {
    std::string treeFileName = outfile;
    treeFileName += ".tre";
    treeStrm = new std::ofstream(treeFileName.c_str(), std::ios::out);
    if (treeStrm == nullptr)
        Msg::error("Could not open tree and/or parameter files");
}

void Mcmc::sample(int n, double lnL) {
    std::string parmFileName = outfile;
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
    std::string treeStr = t->getNewickString();
    *treeStrm << "   tree tree_" << n << " = " << treeStr << "\n";
    
    if (n == numCycles)
        {
        *treeStrm << "end;\n";
        w.closeTSV();
        }
}
