#ifndef MetropolisCoupledMcmc_hpp
#define MetropolisCoupledMcmc_hpp

#include <fstream>
#include <vector>

#include "ThreadPool.hpp"
#include "WriteTSV.hpp"

class PhylogeneticModel;

class MetropolisCoupledMcmc {
    public:
                                            MetropolisCoupledMcmc(void) = delete;
                                            MetropolisCoupledMcmc(unsigned long ng, int pf, int sf, std::vector<PhylogeneticModel*> m);
        void                                run(void);
    
    private:
        //functions
        double                              calcHeating(int idx);
        void                                sample(unsigned long n);
        void                                updateDeltaT(void);
        //objects ordered by memory footprint
        std::deque<bool>                    recentAcceptRej;
        std::vector<PhylogeneticModel*>     models;
        std::vector<double>                 currLnL;
        std::vector<double>                 newLnL;
        std::vector<double>                 currLnP;
        std::vector<double>                 newLnP;
        std::vector<double>                 lnProposalRatio;
        std::vector<double>                 lnLikelihoodRatio;
        std::vector<double>                 lnPriorRatio;
        std::vector<double>                 lnAcceptanceProbabilities;
        std::vector<int>                    indices;
        ThreadPool                          threadPool;
        std::string                         tracerFileName;
        WriteTSV                            w;
        unsigned long                       numCycles;
        double                              deltaT;
        int                                 coldModelIdx;
        int                                 numModels;
        int                                 numSwapsCold;
        int                                 printFrequency;
        int                                 sampleFrequency;
};

#endif
