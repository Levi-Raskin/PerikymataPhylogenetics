#ifndef ParameterIntraspecificVarianceGibbs_hpp
#define ParameterIntraspecificVarianceGibbs_hpp

#include "Eigen/Dense"
#include "Parameter.hpp"
#include "ParameterMatrix.hpp"
#include <string>
#include <vector>

class ParameterIntraspecificMean;

class ParameterIntraspecificVarianceGibbs : public ParameterMatrix {
    public:
                                                        ParameterIntraspecificVarianceGibbs(double prob, std::string n, Eigen::MatrixXd* data, PhylogeneticModel* p);
        double                                          getAcceptanceRatio(void) { return ((double)numAcceptances)/((double)(numAcceptances+numRejections));}
        bool                                            getAdaptiveProposalActive(void);
        const Eigen::MatrixXd&                          getValue(void);
        double                                          lnProbability(void);
        void                                            print(void);
        void                                            setMean(ParameterIntraspecificMean* m) { mean = m; }
        double                                          update(void);
        void                                            updateForAcceptance(void);
        void                                            updateForRejection(void);
        
    private:
        //Objects ordered by memory footprint
        std::vector<Eigen::MatrixXd>                    value;
        Eigen::MatrixXd                                 psi;
        PhylogeneticModel*                              model;
        ParameterIntraspecificMean*                     mean;
        Eigen::MatrixXd*                                tipData;
        double                                          dof;
        double                                          cachedlnP;
        bool                                            useCachedLnP;
        int                                             numAcceptances;
        int                                             numRejections;
        int                                             numtraits;        
        int                                             nObs;        
};

#endif /* ParameterIntraspecificVarianceGibbs_hpp */
