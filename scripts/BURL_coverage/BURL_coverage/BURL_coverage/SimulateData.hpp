#ifndef SimulateData_hpp
#define SimulateData_hpp

#include "Eigen/Dense"
#include "Tree.hpp"
#include <map>

class SimulateData{
    public:
                                                SimulateData(void);
        void                                    checkCredIntdata(void);
        Eigen::MatrixXd*                        getSimulatedData(void){ return &data; };
        std::vector<std::string>                getSimulatedRownames(void){ return rownames; };
        Tree*                                   getSimulatedTree(void){ return tree; };
        double                                  getVCVInCredInt(void);
        double                                  getTipMeanInCredInt(void);
        double                                  getTipVCVInCredInt(void);
        double                                  getImputedInCredInt(void);
        void                                    simulateData(void);
        
    private:
        //coverage instance vars
        Eigen::MatrixXi                         vcvInCredInt;       // ntraits × ntraits
        Eigen::MatrixXi                         tipMeanInCredInt;   // ntips   × ntraits
        std::vector<Eigen::MatrixXi>            tipVCVInCredInt;    // ntips   × ntraits × ntraits
        Eigen::VectorXi                         imputedInCredInt;   // nimp
        
        //prior parameters
        double                                  priorDOF;
        Eigen::MatrixXd                         psi;
        
        //simulated parameters
        Eigen::MatrixXd                         data;
        std::vector<std::string>                rownames;
        Eigen::MatrixXd                         sampledEvoVCV;
        Tree*                                   tree;
        std::map<std::string, Eigen::VectorXd>  trueTipMeans;
        std::map<std::string, Eigen::MatrixXd>  trueTipVCVs;
        std::map<std::string, int>              tipNameToIndex;
        std::vector<std::string>                tipNames;
        std::map<std::pair<int, int>, double>   trueMissingValues;
        
        //misc
        int                                     nimp;
        int                                     nind;
        int                                     nreps;
        int                                     ntips;
        int                                     ntraits;
};

#endif
