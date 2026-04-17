#include "Mcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "WriteTSV.hpp"

#include <iomanip>
#include <iostream>


Mcmc::Mcmc(int ng, int pf, int sf, PhylogeneticModel* m) : numCycles(ng), printFrequency(pf), sampleFrequency(sf), model(m) {
    UserSettings& settings = UserSettings::userSettings();
    tracerFileName = settings.getOutputFile();
    WriteTSV w();
}

void Mcmc::run(void) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double curLnL = model->lnLikelihood();
    double curLnP = model->lnPriorProbability();

    if(numCycles >= std::numeric_limits<unsigned long>::max())
        Msg::error("numCycles requested in greater than largest possible value");

    for (unsigned long n=1; n<=numCycles; n++) {
        double lnProposalRatio = model->update();
        double newLnL = model->lnLikelihood();
        double lnLikelihoodRatio = newLnL - curLnL;
        double newLnP = model->lnPriorProbability();
        double lnPriorRatio = newLnP - curLnP;
        double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

        bool acceptMove = false;
        if (log(rng.uniformRv()) < lnR)
            acceptMove = true;
            
        if (n % printFrequency == 0)
            {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << n << " -- " << curLnL << " -> " << newLnL;
            model->print();
            }
            
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
}

void Mcmc::sample(unsigned long n, double lnL) {
    if(n == 1){
        w.addFilepath(tracerFileName, true);
        std::vector<std::string> cn = {"n", "lnL"};
        std::vector<std::string> headStr = model->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        w.addColumnNamesTSV(cn);
    }

    std::vector<double> dat = {(double)n, lnL};
    std::vector<double> parmStr = model->getParameterString();
    dat.insert( dat.end(), parmStr.begin(), parmStr.end() );
    w.appendDataTSV(dat);

    if (n == numCycles){
        w.closeTSV();
    }
}
