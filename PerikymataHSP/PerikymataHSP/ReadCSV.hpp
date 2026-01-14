#ifndef ReadCSV_hpp
#define ReadCSV_hpp

#include "Eigen/Dense"

#include <string>
#include <vector>

class ReadCSV{
    public:
                                                ReadCSV(std::string filepath, bool rownames, bool colnames);
                                                ReadCSV(std::string filepath, bool rownames, bool colnames, bool string);
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
        Eigen::MatrixXd                         scratch;
        bool                                    readingString;
        std::vector<std::vector<std::string>>   stringData;
};

#endif
