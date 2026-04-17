#ifndef SimulateData_hpp
#define SimulateData_hpp

#include "Eigen/Dense"
#include "Tree.hpp"

class SimulateData{
    public:
                                                SimulateData(void);
        void                                    checkCredInt(void);
        double                                  getVCVInCredInt(void);
        double                                  getTipMeanInCredInt(void);
        double                                  getTipVCVInCredInt(void);
        double                                  getImputedInCredInt(void);
        void                                    simulateData(void);
        
    private:
        Eigen::MatrixXi                         vcvInCredInt;       // ntraits × ntraits
        Eigen::MatrixXi                         tipMeanInCredInt;   // ntips   × ntraits
        std::vector<Eigen::MatrixXi>            tipVCVInCredInt;    // ntips   × ntraits × ntraits
        Eigen::MatrixXi                         imputedInCredInt;   // nimp
        
        //prior parameters
        double                                  priorDOF;
        Eigen::MatrixXd                         psi;
        
        //simulated parameters
        std::vector<Eigen::MatrixXd>            sampledTipVCV;
        Eigen::MatrixXd                         sampledEvoVCV;
        Tree                                    tree;
        std::map<std::string, Eigen::VectorXd>  trueTipMeans;
        std::map<std::string, Eigen::MatrixXd>  trueTipVCVs;
        std::map<std::string, int>              tipNameToIndex;
        Eigen::VectorXd                         trueMissingValues;
};

#endif
