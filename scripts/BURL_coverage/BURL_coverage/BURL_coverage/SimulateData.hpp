#ifndef SimulateData_hpp
#define SimulateData_hpp
#include "Eigen/Dense"
#include "Tree.hpp"
#include <chrono>
#include <fstream>
#include <map>
class SimulateData{
    public:
                                                SimulateData(void);
        void                                    checkCredInt(void);
        void                                    endIncrement(void);
        Eigen::MatrixXd*                        getSimulatedData(void){ return &data; };
        std::vector<std::string>                getSimulatedRownames(void){ return rownames; };
        Tree*                                   getSimulatedTree(void){ return tree; };
        double                                  getIncrementedElapsed(void) const { return incrementedElapsed; };
        void                                    print(void);
        void                                    startIncrement(void);
        void                                    writeCoverage(void);
        void                                    simulateData(void);
        
    private:
        std::pair<double,double>                meanAndSE(const std::vector<double>& v) const;
        //rank-uniformity helpers
        void                                    checkRankUniformity(const Eigen::MatrixXd& rDat, const std::vector<std::string>& cn);
        Eigen::VectorXd                         postBurninThinned(const Eigen::VectorXd& col) const;
        Eigen::VectorXd                         vcvTraceTrace(const Eigen::MatrixXd& rDat, const std::vector<std::string>& cn,
                                                                const std::string& vcvPrefix, int nTraits) const;
        int                                     computeRank(const Eigen::VectorXd& posterior, double trueVal) const;
        void                                    writeRankRow(const std::string& label, int rank, int L);
        
        //coverage instance vars
        Eigen::MatrixXi                         vcvInCredInt;       // ntraits × ntraits
        Eigen::MatrixXi                         tipMeanInCredInt;   // ntips   × ntraits
        std::vector<Eigen::MatrixXi>            tipVCVInCredInt;    // ntips   × ntraits × ntraits
        Eigen::VectorXi                         imputedInCredInt;   // nimp
        
        std::vector<double>                     evoVCVRepCoverage;
        std::vector<double>                     tipVCVRepCoverage;
        std::vector<double>                     tipMeanRepCoverage;
        std::vector<double>                     missingRepCoverage;
        std::vector<double>                     totalRepCoverage;
        
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
        
        //timing: accumulates wall-clock time ONLY while an MCMC/MC3 run is
        //active, i.e. only between startIncrement()/endIncrement() pairs --
        //excludes simulateData(), checkCredInt(), and print() overhead.
        std::chrono::steady_clock::time_point   incrementStart;
        double                                  incrementedElapsed;
        bool                                    incrementRunning;
        
        //rank-uniformity output
        std::ofstream                           rankOut;
        bool                                    rankFileInitialized;
        double                                  rankBurninFraction;   // fraction of trace ROWS discarded before ranking
        int                                     rankThinStride;       // keep every Nth retained row (1 = no extra thinning)
        
        //misc
        int                                     nimp;
        int                                     nind;
        int                                     nreps;
        int                                     ntips;
        int                                     ntraits;
        int                                     trials;
};
#endif
