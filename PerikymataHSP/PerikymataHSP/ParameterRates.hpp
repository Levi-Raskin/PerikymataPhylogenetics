#ifndef ParameterRates_hpp
#define ParameterRates_hpp

#include "Eigen/Dense"
#include "Parameter.hpp"
#include <deque>
#include <string>
#include <vector>

class PhylogeneticModel;

class ParameterRates : public Parameter {
    public:
                                        ParameterRates(double prob, PhylogeneticModel* p, std::string n, int nt);
        double                          getAcceptanceRatio(void) { return ((double)numAcceptances)/((double)(numAcceptances+numRejections));}
        bool                            getAdaptiveProposalActive(void) { return adaptiveProposalActive; }
        Eigen::VectorXd                 getValue(void) { return value[0].array().exp(); } // 0 is the one we update, 1 is the one we don't (last currently accepted value); on log scale
        double                          lnProbability(void);
        void                            print(void);
        void                            setValue(Eigen::VectorXd x) { value[0] = value[1] = x.array().log(); }
        double                          update(void);
        void                            updateForAcceptance(void);
        void                            updateForRejection(void);
    private:
        double                          updateMH(void);
        double                          updateHMC(void);
        Eigen::VectorXd                 a; //dirichlet param
        int                             nTraits;
        int                             numAcceptances;
        int                             numRejections;
        std::vector<Eigen::VectorXd>    value;
        double                          windowSize;
        std::deque<bool>                recentAcceptRej;

};

#endif /* ParameterRates_hpp */
