#ifndef Probability_hpp
#define Probability_hpp

#ifndef PI
#    define PI 3.141592653589793
#endif

#include "Eigen/Dense"
#include <vector>
class RandomVariable;
////eigen matrixXd forward dec
//namespace Eigen {
//    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    class Matrix;
//    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
//}

thread_local static bool availableNormalRv = false;
thread_local static double extraNormalRv = 0.0;



namespace  Probability {

    namespace Beta {
    
        double  pdf(double alpha, double beta, double x);
        double  lnPdf(double alpha, double beta, double x);
        double  rv(RandomVariable* rng, double a, double b);
        double  cdf(double alpha, double beta, double x);
        double  quantile(double alpha, double beta, double x);
    }

    namespace ChiSquare {
    
        double  pdf(double v, double x);
        double  lnPdf(double v, double x);
        double  rv(RandomVariable* rng, double v);
        double  cdf(double v, double x);
        double  quantile(double prob, double v);
    }

    namespace Dirichlet {
    
        double  pdf(const std::vector<double> &a, const std::vector<double> &z);
        double  lnPdf(const std::vector<double> &a, const std::vector<double> &z);
        double  lnPdf(const Eigen::VectorXd &a, const Eigen::VectorXd &z);
        void    rv(RandomVariable* rng, const std::vector<double> &a, std::vector<double> &z);
        void    rv(RandomVariable* rng, const Eigen::VectorXd &a, Eigen::VectorXd &z);
    }

    namespace Gamma {
    
        double  pdf(double alpha, double beta, double x);
        double  lnPdf(double alpha, double beta, double x);
        double  rv(RandomVariable* rng, double alpha, double beta);
        double  cdf(double alpha, double beta, double x);
        void    discretization(std::vector<double> &catRate, double a, double b, int nCats, bool median);
    }
    
    namespace Geometric {
    
        int     rv(RandomVariable* rng, double p);
    }

    namespace Exponential {
    
        double  pdf(double lambda, double x);
        double  lnPdf(double lambda, double x);
        double  rv(RandomVariable* rng, double lambda);
        double  cdf(double lambda, double x);
    }

    namespace Normal {
    
        double  pdf(double mu, double sigma, double x);
        double  lnPdf(double mu, double sigma2, double x);
        double  rv(RandomVariable* rng) ;
        double  rv(RandomVariable* rng, double mu, double sigma) ;
        double  cdf(double mu, double sigma, double x);
        double  invCdf(double mu, double sigma, double scope);
        double  invCdfStandardNormal(double scope);
        double  quantile(double mu, double sigma, double p);
    }
    
    namespace TruncatedNormal {
        
        double  lnPdf(double x, double mu, double sigma, double lower, double upper);
        double  rv(RandomVariable* rng, double mu, double sigma, double lower, double upper) ;
    }

    namespace Uniform {
    
        double  pdf(double low, double high, double x);
        double  lnPdf(double low, double high, double x);
        double  rv(RandomVariable* rng, double low, double high);
        double  cdf(double low, double high, double x);
    }
    
    namespace MultivariateT{
        double          lnPdf(Eigen::MatrixXd* dat, Eigen::MatrixXd* u, Eigen::MatrixXd* scale, double dof);
    }

    namespace MultivariateNormal{
        std::vector<double> rv(RandomVariable* rng, std::vector<double> u, Eigen::MatrixXd* VCV);
        Eigen::VectorXd rv(RandomVariable* rng, Eigen::VectorXd u, Eigen::MatrixXd* VCV);
        Eigen::VectorXd rv(RandomVariable* rng, const Eigen::VectorXd& u, const Eigen::MatrixXd& VCV);
        void            rv(RandomVariable* rng,
                            Eigen::VectorXd& out,
                            const Eigen::VectorXd& means,
                            const Eigen::MatrixXd& vcvCholLower);
        double          lnPdf(Eigen::MatrixXd* vec, Eigen::MatrixXd* distribVec, Eigen::MatrixXd* VCV);
        double          lnPdf(Eigen::VectorXd* vec, Eigen::VectorXd* u, Eigen::MatrixXd* VCV);
        double          lnPrecisionMatrixPdf(Eigen::VectorXd* vec, Eigen::VectorXd* u, Eigen::MatrixXd* VCVInv);
        std::pair<Eigen::MatrixXd, Eigen::MatrixXd> productMVN(Eigen::MatrixXd& mean0, Eigen::MatrixXd& var0,Eigen::MatrixXd& mean1, Eigen::MatrixXd& var1);
    }
    
    namespace Wishart{
        Eigen::MatrixXd rv(RandomVariable* rng, Eigen::MatrixXd* V, double n);
        double          lnPdf(Eigen::MatrixXd* support, Eigen::MatrixXd* scale, double dof);
    }
    
    namespace InverseWishart{
        Eigen::MatrixXd rv(RandomVariable* rng, const Eigen::MatrixXd& psi, double n);
        void            rv(RandomVariable* rng,
                            Eigen::MatrixXd& out,
                            const Eigen::MatrixXd& psiInvCholLower,
                            double n);            
        Eigen::MatrixXd rv(RandomVariable* rng, Eigen::MatrixXd* psi, double n);
        double          lnPdf(const Eigen::MatrixXd& x, const Eigen::MatrixXd& psi, double v);
        double          lnPdf(Eigen::MatrixXd* support, Eigen::MatrixXd* scale, double dof);
    }
    
    namespace LKJ{
        double          lnPdf(Eigen::MatrixXd* support, double shape);
    }
    
    namespace Helper {
        
        double  beta(double a, double b);
        double  chebyshevEval(double x, const double *a, const int n);
        double  epsilon(void);
        double  gamma(double x);
        bool    isFinite(double x);
        double  lnFactorial(int n);
        double  lnBeta(double a, double b);
        double  lnGamma(double a);
        double  lnGammacor(double x);
        double  incompleteBeta(double a, double b, double x);
        double  incompleteGamma(double x, double alpha, double LnGamma_alpha);
        void    normalize(std::vector<double>& vec);
        double  pointNormal(double prob);
        double  rndGamma(RandomVariable* rng, double s);
        double  rndGamma1(RandomVariable* rng, double s);
        double  rndGamma2(RandomVariable* rng, double s);
    }
}

#endif
