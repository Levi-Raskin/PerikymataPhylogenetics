#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <fstream>
class PhylogeneticModel;

class Mcmc {
    public:
                                Mcmc(void) = delete;
                                Mcmc(int ng, int pf, int sf, PhylogeneticModel* m);
                               ~Mcmc(void);
        void                    run(void);
        void                    run(int nc);
    
    private:
        void                    openFiles(void);
        void                    sample(unsigned long n, double lnL);
        void                    sample(unsigned long n, double lnL, double lnP);
        int                     numCycles;
        int                     printFrequency;
        int                     sampleFrequency;
        PhylogeneticModel*      model;
        std::ofstream*          treeStrm;
};

#endif
