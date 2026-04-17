#ifndef ParameterIntraspecificMean_hpp
#define ParameterIntraspecificMean_hpp

#include "Eigen/Dense"
#include "Parameter.hpp"
#include "ParameterMatrix.hpp"
#include <string>
#include <vector>

class ParameterIntraspecificVariance;
class PhylogeneticModel;

class ParameterIntraspecificMean : public Parameter {
    public:
                                                        ParameterIntraspecificMean(double prob, std::string n, Eigen::MatrixXd* data, PhylogeneticModel* p);
        double                                          getAcceptanceRatio(void){
                                                            const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
                                                            const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
                                                            return acceptanceRate;
                                                        }
        bool                                            getAdaptiveProposalActive(void) { return adaptiveProposalActive; }
        const Eigen::VectorXd&                          getValue(void) { return mean[0]; }
        double                                          lnProbability(void);
        void                                            print(void);
        void                                            resetAdaptiveTuning(void) {
                                                            numAcceptances = 0;
                                                            numRejections = 0;
                                                            recentAcceptRej.clear();
                                                            ngAdaptive = 50000;
                                                            windowSize = 0.1;
                                                        }
        void                                            setNGAdaptive(int n) { ngAdaptive = n; }
        void                                            setVarianceCovarianceMatrix(ParameterMatrix* vcv) { varianceCovariance = vcv; }
        double                                          update(void);
        void                                            updateForAcceptance(void);
        void                                            updateForRejection(void);
    private:
        //functions
        double                                          updateMHmvNDraw(void);
        //objects oriented in decreasing memory impact
        std::vector<Eigen::VectorXd>                    mean;
        Eigen::MatrixXd                                 psi;
        Eigen::MatrixXd                                 proposalCov;
        Eigen::MatrixXd                                 proposalCholLower;
        Eigen::VectorXd                                 mu0;
        std::deque<bool>                                recentAcceptRej;
        PhylogeneticModel*                              model;
        Eigen::MatrixXd*                                tipData;
        ParameterMatrix*                                varianceCovariance;
        double                                          windowSize;
        double                                          stepSize;
        double                                          dof;
        double                                          lambda;
        double                                          targetAcceptanceRate;
        double                                          upperAcceptanceRate;
        double                                          lowerAcceptanceRate;
        int                                             ngAdaptive;
        int                                             nObs;
        int                                             nTraits;
        int                                             numAcceptances;
        int                                             numRejections;
        //rearrange
        Eigen::MatrixXd                                 empiricalCovariance;
        Eigen::VectorXd                                 empiricalMean;
        std::deque<Eigen::VectorXd>                     recentSamples;
        int                                             covarianceUpdateFreq;
        bool                                            useEmpiricalCovariance;
        void                                            updateEmpiricalCovariance(void);
};

#endif
