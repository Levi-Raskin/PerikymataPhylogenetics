#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <fstream>
class LandmarkAlignmentModel;
class MultivariateNormalPhylogeneticModel;
class PhylogeneticModel;

class Mcmc {
    public:
                                Mcmc(void) = delete;
                                Mcmc(int ng, int pf, int sf, PhylogeneticModel* m, std::string of);
                               ~Mcmc(void);
        void                    run(void);
    
    private:
        void                    openFiles(void);
        void                    sample(int n, double lnL);
        int                     nChains; //assuming parallized
        int                     numCycles;
        std::string             outfile;
        int                     printFrequency;
        int                     sampleFrequency;
        PhylogeneticModel*      model;
        std::ofstream*          treeStrm;
};

#endif
