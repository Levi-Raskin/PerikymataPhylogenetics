#ifndef Mcmc_hpp
#define Mcmc_hpp

#include "WriteTSV.hpp"
class PhylogeneticModel;

class Mcmc {
    public:
                                Mcmc(void) = delete;
                                Mcmc(int ng, int pf, int sf, PhylogeneticModel* m);
        void                    run(void);
    
    private:
        void                    sample(unsigned long n, double lnL);
        PhylogeneticModel*      model;
        WriteTSV                w;
        std::string             tracerFileName;
        int                     numCycles;
        int                     printFrequency;
        int                     sampleFrequency;
};

#endif
