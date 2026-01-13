#ifndef ParameterInvVCV_hpp
#define ParameterInvVCV_hpp

#include "Eigen/Dense"
#include "ParameterMatrix.hpp"
#include <deque>
#include <string>

class ParameterCorrelationMatrix;
class ParameterDouble;
class ParameterRates;
class PhylogeneticModel;

class ParameterInvVCV : public ParameterMatrix{
    public:
                                    ParameterInvVCV(double prob, std::string n, int numberOfTraits, PhylogeneticModel* p);
        double                      getAcceptanceRatio(void);
        bool                        getAdaptiveProposalActive(void);
        Eigen::MatrixXd             getValue(void); //returns VCV-1
        double                      lnProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        
    private:
        std::vector<Parameter*>     components;
        ParameterCorrelationMatrix* correlationMatrix;
        int                         numtraits;
        ParameterRates*             rates;
        double                      targetAcceptanceRate;
        Parameter*                  updatedComponent;
};

#endif /* ParameterInvVarianceCovarianceMatrixRBDecomp_hpp */
