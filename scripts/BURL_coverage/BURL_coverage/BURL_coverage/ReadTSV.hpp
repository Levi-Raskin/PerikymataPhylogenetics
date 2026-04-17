#ifndef ReadTSV_hpp
#define ReadTSV_hpp

#include <string>
#include <vector>

namespace Eigen {
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    class Matrix;
    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
}

class ReadTSV{
    public:
                                                ReadTSV(std::string filepath, bool rownames, bool colnames);
                                                ReadTSV(std::string filepath, bool rownames, bool colnames, bool string);
        std::vector<std::string>                getColnames(void){return colnames;}
        Eigen::MatrixXd                         getEigenMat(void);
        std::vector<std::vector<double>>        getReadData(void){return doubleData;}
        std::vector<std::vector<std::string>>   getReadStringData(void){return stringData;}
        std::vector<std::string>                getRownames(void){return rownames;}
        void                                    print(void);
    private:
        std::vector<std::vector<double>>        doubleData;
        std::vector<std::string>                colnames;
        std::vector<std::string>                rownames;
        bool                                    hasColnames;
        bool                                    hasRownames;
        bool                                    readingString;
        std::vector<std::vector<std::string>>   stringData;
};

#endif 
