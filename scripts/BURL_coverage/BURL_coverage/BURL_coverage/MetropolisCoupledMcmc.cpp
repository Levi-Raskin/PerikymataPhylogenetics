#include "MetropolisCoupledMcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "ThreadPool.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "WriteTSV.hpp"

#include <future>
#include <iomanip>
#include <iostream>

MetropolisCoupledMcmc::MetropolisCoupledMcmc(unsigned long ng, int pf, int sf, std::vector<PhylogeneticModel*> m) : numCycles(ng), printFrequency(pf),
    sampleFrequency(sf),
    models(m),
    numModels(models.size()),
    coldModelIdx(-1),
    numSwapsCold(0),
    deltaT(0.2),
    threadPool(models.size()){
    currLnL.reserve(numModels); //needs to be reserve for pushback
    newLnL.resize(numModels);
    currLnP.reserve(numModels); //needs to be reserve for pushback
    newLnP.resize(numModels);
    indices.reserve(numModels); //needs to be reserve for pushback
    lnProposalRatio.resize(numModels);
    lnLikelihoodRatio.resize(numModels);
    lnPriorRatio.resize(numModels);
    lnAcceptanceProbabilities.resize(numModels);
    
    UserSettings& settings = UserSettings::userSettings();
    tracerFileName = settings.getOutputFile();
    WriteTSV w();
}

double MetropolisCoupledMcmc::calcHeating(int idx){
    if(idx == 0)
        return 1.0;
    return (1 / (1 + deltaT * idx));
}

void MetropolisCoupledMcmc::run(void) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    int idx = 0;
    for(PhylogeneticModel* m : models){
        currLnL.push_back(m->lnLikelihood());
        currLnP.push_back(m->lnPriorProbability());
        indices.push_back(idx);
        idx++;
    }
    
    if(numCycles >= std::numeric_limits<unsigned long>::max())
        Msg::error("Numer of cycles requested is greater than largest possible value");

    for (unsigned long n=1; n<=numCycles; n++){
        if(n < 10000 && n % 50 == 0)
            updateDeltaT();
        
        std::vector<std::future<void>> futures;
        futures.reserve(numModels);

        for(int i = 0; i < numModels; i++){
            auto promise = std::make_shared<std::promise<void>>();
            futures.push_back(promise->get_future());

            threadPool.enqueue([this, i, promise]() {
                lnProposalRatio[i] = models[i]->update();
                newLnL[i]          = models[i]->lnLikelihood();
                newLnP[i]          = models[i]->lnPriorProbability();

                lnLikelihoodRatio[i] = newLnL[i] - currLnL[i];
                lnPriorRatio[i]      = newLnP[i] - currLnP[i];

                int    chainIdx = indices[i];
                double heat     = (chainIdx == 0) ? 1.0 : calcHeating(chainIdx);
                lnAcceptanceProbabilities[i] = heat * (lnLikelihoodRatio[i] + lnPriorRatio[i])
                                               + lnProposalRatio[i];

                promise->set_value();
            });
        }

        // Wait for all chain updates to complete before proceeding
        for(auto& f : futures)
            f.get();


        // accept or reject the proposed state
        coldModelIdx = -1;
        for (int i = 0; i < numModels; i++) {
            if (indices[i] == 0) coldModelIdx = i;
            if (log(rng.uniformRv()) < lnAcceptanceProbabilities[i]) {
                currLnL[i] = newLnL[i];
                currLnP[i] = newLnP[i];
                models[i]->updateForAcceptance();
            } else {
                models[i]->updateForRejection();
            }
        }
        
        if(coldModelIdx == -1)
            Msg::error("did not find cold model");
        
        if (n % printFrequency == 0){
            std::cout << std::fixed << std::setprecision(2);
            const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
            const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
            std::cout << n << " -- "
                << currLnL[coldModelIdx] << " -> " << newLnL[coldModelIdx]
                << " | Chain swap acceptance ratio: " << acceptanceRate << "\n";
            models[coldModelIdx]->print();
        }
            
        //choose two chains and swap
        int chain1 = (int)(rng.uniformRv() * numModels);
        int chain2 = -1;
        do{
            chain2 = (int)(rng.uniformRv() * numModels);
        }while(chain1 == chain2);
    
        int idx1 = indices[chain1];
        int idx2 = indices[chain2];
        
        double chain1LnPost = currLnL[chain1] + currLnP[chain1];
        double chain2LnPost = currLnL[chain2] + currLnP[chain2];
        
        double chain1Heating = calcHeating(idx1);
        double chain2Heating = calcHeating(idx2);
        
        double lnProbSwap = chain2Heating * chain1LnPost + chain1Heating * chain2LnPost;
        double lnProbStay = chain1Heating * chain1LnPost + chain2Heating * chain2LnPost;
        double lnAcceptanceSwap =lnProbSwap - lnProbStay;
        
        //accept or reject swap
        if(log(rng.uniformRv()) < lnAcceptanceSwap){
            indices[chain1] = idx2;
            indices[chain2] = idx1;
            recentAcceptRej.push_back(true);
            if(idx2 == 0 || idx1 == 0)
                numSwapsCold++;
        }else{
            recentAcceptRej.push_back(false);
        }
        
        if (recentAcceptRej.size() > 10000)
            recentAcceptRej.pop_front();
        
        if (n == 1 || n == numCycles || n % sampleFrequency == 0 )
            sample(n);
    }
}

void MetropolisCoupledMcmc::sample(unsigned long n) {
    if(n == 1){
        w.addFilepath(tracerFileName, true);
        std::vector<std::string> cn = {"n", "lnL"};
        std::vector<std::string> headStr = models[coldModelIdx]->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        w.addColumnNamesTSV(cn);
    }
        
    std::vector<double> dat = {(double)n, currLnL[coldModelIdx]};
    std::vector<double> parmStr = models[coldModelIdx]->getParameterString();
    dat.insert( dat.end(), parmStr.begin(), parmStr.end() );
    w.appendDataTSV(dat);
    
    if (n == numCycles){
        w.closeTSV();
    }
}

void MetropolisCoupledMcmc::updateDeltaT(void) {
    if(recentAcceptRej.size() > 50){
        const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
        const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
        
        constexpr double targetAcceptance = 0.23;
        constexpr double lowerAcceptance = targetAcceptance - 0.1;
        constexpr double upperAcceptance = targetAcceptance + 0.1;
        if(acceptanceRate < lowerAcceptance){
            deltaT *= 0.99;
        }else if (acceptanceRate > upperAcceptance){
            deltaT *= 1.01;
        }
    }
}
