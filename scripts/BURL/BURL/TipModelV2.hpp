#ifndef TipModelV2_hpp
#define TipModelV2_hpp

#include "PhylogeneticModel.hpp"
#include "Eigen/Dense"

class ParameterIntraspecificMean;
class ParameterIntraspecificVariance;
class ParameterIntraspecificVarianceGibbs;
class ParameterDouble;

class TipModelV2 : public PhylogeneticModel{
    public:
                                                                            TipModelV2(std::string tn, Eigen::MatrixXd d, PerikymataHSPv4* m);
                                                                           ~TipModelV2(void);
        const Eigen::VectorXd&                                              getTipMean(void);
        std::vector<std::string>                                            getParameterNames(void);
        std::vector<double>                                                 getParameterString(void);
        double                                                              lnLikelihood(void);
        double                                                              lnPriorProbability(void);
        void                                                                print(void);
        double                                                              update(void);
        void                                                                updateForAcceptance(void);
        void                                                                updateForRejection(void);
        
    private:
        double                                                              computeLnLikelihood(void);
        double                                                              computeLnPriorProbability(void);
        void                                                                updateTipDataComplete(void);
        void                                                                updatePkGibbs(void);
        Eigen::LLT<Eigen::MatrixXd>                                         sigmaChol;
        Eigen::LLT<Eigen::MatrixXd>                                         sigma22Solver;
        Eigen::MatrixXd                                                     tipDataIncomplete; //raw data; stays pristene
        Eigen::MatrixXd                                                     tipDataComplete; //used to calc. likelihood of intraspecific parms
        Eigen::MatrixXd                                                     sigma11;
        Eigen::MatrixXd                                                     sigma12;
        Eigen::MatrixXd                                                     sigma21;
        Eigen::MatrixXd                                                     sigma22;
        Eigen::MatrixXd                                                     sigma22Inv_sigma21;
        Eigen::MatrixXd                                                     sigmaCond;
        Eigen::MatrixXd                                                     sigma;
        Eigen::MatrixXd                                                     L;
        Eigen::MatrixXd                                                     scratchMat;
        Eigen::VectorXd                                                     u1;
        Eigen::VectorXd                                                     u2;
        Eigen::VectorXd                                                     x2;
        Eigen::VectorXd                                                     xDiff;
        Eigen::VectorXd                                                     x2_minus_u2;
        Eigen::VectorXd                                                     sigma22Inv_diff;
        Eigen::VectorXd                                                     uCond;
        Eigen::VectorXd                                                     newVals;
        Eigen::VectorXd                                                     mu;
        Eigen::VectorXd                                                     scratchVec;
        std::map<std::pair<int, int>, ParameterDouble*>                     missingPkVals;
        std::vector<std::string>                                            parmNames;
        std::vector<double>                                                 parmValues;
        std::vector<int>                                                    missingIndices;
        std::vector<int>                                                    obsIndices;
        std::vector<ParameterDouble*>                                       updatedImpPkDoubles;
        PhylogeneticModel*                                                  model;
        ParameterIntraspecificMean*                                         taxonMean;
        ParameterIntraspecificVarianceGibbs*                                taxonVariance;
        std::string                                                         tipName;
        double                                                              sigmaLogDet;
        double                                                              cachedLnL;
        double                                                              cachedLnP;
        double                                                              log2pi;
        double                                                              numCnumR;
        double                                                              term1;
        bool                                                                lnLDirty;
        bool                                                                lnPDirty;
        bool                                                                hasMissingData;
        bool                                                                gibbsPkUpdate;
        int                                                                 numCols;
        int                                                                 numRows;
        int                                                                 numImputationRejections;
        int                                                                 numImputationAcceptances;
};

#endif /* TipModelV2_hpp */
