#ifndef Utility_hpp
#define Utility_hpp

#include <utility>
#include <vector>

//eigen matrixxd forward declaration
namespace Eigen {
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    class Matrix;
    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
    using VectorXd = Matrix<double, -1, 1, 0, -1, 1>;    // Dynamic size double column vector
}

namespace  Utility {
    namespace Bayesian{
        std::pair<double, double>   credibleInterval(Eigen::VectorXd v);
        std::pair<double, double>   credibleIntervalBurnIn(Eigen::VectorXd v, double bi);
        std::pair<double, double>   hpdIntervalBurnIn(Eigen::VectorXd v, double bi);
        bool                        multimodalCredibleInterval(double x, Eigen::VectorXd v, double bi);
    }
    namespace Shapes{
        Eigen::MatrixXd generateUnitCirclePoints(int n);
        Eigen::MatrixXd generateUnitSpherePoints(int n);
    }
    namespace EigenUtils {
        void            deleteColumn(Eigen::MatrixXd& matrix, unsigned int colToDelete);
        void            deleteRow(Eigen::MatrixXd& matrix, unsigned int rowToDelete);
        Eigen::MatrixXd kroneckerProduct(Eigen::MatrixXd& matrix0, Eigen::MatrixXd& matrix1);
        void            printEigen(Eigen::MatrixXd e);
        void            printEigenR(Eigen::MatrixXd e);
        Eigen::MatrixXd vectorMatrix2Eigen(std::vector<std::vector<double>> v);
    }
}

#endif
