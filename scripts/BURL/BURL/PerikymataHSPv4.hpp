#ifndef PerikymataHSPv4_hpp
#define PerikymataHSPv4_hpp

#include "Eigen/Dense"
#include "MultivariateBrownianMotionV2.hpp"
#include <string>

class ParameterDouble;
class ParameterMatrix;
class ParameterRates;
class ParameterVectorDouble;
class TipModelV2;
class TipModel;
class Tree;

class PerikymataHSPv4 : public MultivariateBrownianMotionV2{
    public:
                                                                                    PerikymataHSPv4(Tree* backbone, std::vector<std::string> datRN, Eigen::MatrixXd* dat);
                                                                                   ~PerikymataHSPv4(void);
        std::vector<std::string>                                                    getParameterNames(void) override;
        std::vector<double>                                                         getParameterString(void) override;
        double                                                                      lnLikelihood(void) override;
        double                                                                      lnPriorProbability(void) override;
        void                                                                        print(void) override;
        double                                                                      update(void) override;
        void                                                                        updateForAcceptance(void) override;
        void                                                                        updateForRejection(void) override;
    private:
        //functions
        void                                                                        calculateTipMeans(void);
        void                                                                        updatePkGibbs(void);
        //objects in order of memory footprint
        Eigen::MatrixXd                                                             scratch;
        Eigen::MatrixXd                                                             tipMeansConcat;
        std::unordered_map<std::string, TipModelV2*>                                tipModels;
        std::unordered_map<std::string, int>                                        tipIdxs;
        TipModelV2*                                                                 updatedTipModel;
        std::vector<std::string>                                                    tipNames;
        std::vector<double>                                                         parmValues;
        std::vector<double>                                                         scratchVec;
        Tree                                                                        fixedTree;
        bool                                                                        tipUpdate;
        bool                                                                        updateTipsOn;
};

#endif
