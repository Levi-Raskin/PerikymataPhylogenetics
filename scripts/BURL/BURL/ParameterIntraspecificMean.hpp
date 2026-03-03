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
                                                            rwScaleFactor = 2.38 / std::sqrt((double)nTraits * (double)nObs);
                                                            malaEpsilon = 1.0 / (std::pow((double)nTraits, 1.0/3.0) * std::sqrt((double)nObs));
                                                            malaDriftCoeff = 0.5 * malaEpsilon * malaEpsilon * (double)nObs;
                                                        }
        void                                            setNGAdaptive(int n) { ngAdaptive = n; }
        void                                            setUsePosteriorDraw(bool b) { usePosteriorDraw = b; }
        void                                            setVarianceCovarianceMatrix(ParameterMatrix* vcv) { varianceCovariance = vcv; }
        double                                          update(void);
        double                                          updatePosteriorDraw(void);
        void                                            updateForAcceptance(void);
        void                                            updateForRejection(void);
    private:
        //functions
        double                                          updateMHSingleElement(void);
        double                                          updateMHmvNDraw(void);
        double                                          updateConditionalDraw(void);
        double                                          updateHMC(void);
        double                                          updateHMCSingleElement(void);
        double                                          updatePreconditionedRW(void);
        double                                          updateDataLikelihoodMALA(void);
        double                                          calculatePosterior(void);
        double                                          numericalGradientSingleElement(int idx);
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
        bool                                            usePosteriorDraw;
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
        // Preconditioned RW proposal members
        Eigen::MatrixXd                                 cholSigma;
        Eigen::VectorXd                                 zDraw;
        double                                          rwScaleFactor;
        // Data-likelihood MALA proposal members
        Eigen::VectorXd                                 dataMean;           // x̄ᵢ cached
        Eigen::VectorXd                                 driftOld;           // h(µ_old)
        Eigen::VectorXd                                 driftNew;           // h(µ_new)
        Eigen::VectorXd                                 residFwd;           // µ_new - µ_old - h(µ_old)
        Eigen::VectorXd                                 residRev;           // µ_old - µ_new - h(µ_new)
        double                                          malaEpsilon;        // MALA step size ε
        double                                          malaDriftCoeff;     // (ε²/2)·nᵢ precomputed


};

#endif
