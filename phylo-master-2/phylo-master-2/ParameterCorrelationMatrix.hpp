#ifndef ParameterCorrelationMatrix_hpp
#define ParameterCorrelationMatrix_hpp

#include "Eigen/Dense"
#include "Parameter.hpp"
#include <deque>
#include <string>
#include <vector>

class ParameterCorrelationMatrix : public Parameter {
    public:
                                        ParameterCorrelationMatrix(double prob, PhylogeneticModel* p, std::string n, int nt);
        double                          getAcceptanceRatio(void) { return ((double)numAcceptances)/((double)(numAcceptances+numRejections));}
        bool                            getAdaptiveProposalActive(void) { return adaptiveProposalActive; }
        Eigen::MatrixXd                 getValue(void) {return value[0];}// 0 is the one we update, 1 is the one we don't (last currently accepted value
        double                          lnProbability(void);
        void                            print(void);
        double                          update(void);
        void                            updateForAcceptance(void);
        void                            updateForRejection(void);
    private:
        int                             nTraits;
        int                             numAcceptances;
        int                             numRejections;
        std::vector<Eigen::MatrixXd>    value;
        double                          windowSize;
        std::deque<bool>                recentAcceptRej;
};

#endif
