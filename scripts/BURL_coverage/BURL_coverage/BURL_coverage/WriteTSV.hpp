#ifndef WriteTSV_hpp
#define WriteTSV_hpp

#include "Eigen/Dense"

#include <fstream>
#include <string>
#include <vector>

class WriteTSV{
    public:
                                    WriteTSV(void);
                                    WriteTSV(std::string filepath, bool overwrite);
                                   ~WriteTSV(void);
        void                        addColumnNamesTSV(std::vector<std::string> cn);
        void                        addRownamesTSV(std::vector<std::string> rn);
        void                        addFilepath(std::string fp, bool overwrite);
        void                        appendDataTSV(Eigen::MatrixXd x);
        void                        appendDataTSV(std::vector<double> data);
        void                        appendDataTSV(std::vector<std::string> data);
        void                        appendDataTSV(std::string data);
        void                        closeTSV(void);
        void                        setRNFlag(bool b) { rnFlag = b;}
        
    private:
        std::fstream                fout;
        int                         numCols;
        int                         numRows;
        std::string                 filepath;
        std::vector<std::string>    rownames;
        bool                        rnFlag;
};

#endif
