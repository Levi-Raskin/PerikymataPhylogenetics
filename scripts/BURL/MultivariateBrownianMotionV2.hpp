//
//  MultivariateBrownianMotionV2.hpp
//  phylo-master
//
//  Created by Levi Raskin on 2/28/26.
//

#ifndef MultivariateBrownianMotionV2_hpp
#define MultivariateBrownianMotionV2_hpp


#include "Eigen/Dense"
#include "PhylogeneticModel.hpp"
#include "ParameterVarianceCovarianceMatrixV3.hpp"
#include "Tree.hpp"
#include <map>
#include <string>
#include <vector>

class Node;
class Parameter;
class ParameterDouble;
class ParameterMatrix;
class ParameterTree;

class MultivariateBrownianMotionV2 : public PhylogeneticModel{
    public:
                                                                            MultivariateBrownianMotionV2();
                                                                           ~MultivariateBrownianMotionV2(void);
        void                                                                addData(std::vector<std::string> rn, Eigen::MatrixXd* data);
        std::vector<std::string>                                            getParameterNames(void);
        std::vector<double>                                                 getParameterString(void);
        double                                                              lnLikelihood(void);
        double                                                              lnPriorProbability(void);
        void                                                                print(void);
        void                                                                setTree(Tree* t);
        double                                                              update(void);
        void                                                                updateForAcceptance(void);
        void                                                                updateForRejection(void);

    protected:
        //functions
        double                                                              calculatePosteriorProbability(void);
        void                                                                instantiateIndependentContrasts(void);
        void                                                                normalizeProposalProbability(void);
        void                                                                populateObservedData();
        void                                                                updateVarianceCovarianceMatrix(void);
        
        // Tip data
        // indexed rowise by tip index-- tip 0 corresponds to row 0 etc.
        std::vector<Eigen::VectorXd>                                        nodeVals;
        std::vector<Eigen::VectorXd>                                        contrasts;
        std::vector<std::pair<int, double>>                                 modifiedBranches;
        std::vector<double>                                                 branchLength;
        std::vector<std::string>                                            originalDataRownames;
        std::vector<Node*>                                                  dpseq;
        std::vector<Node*>                                                  nDesc;
        Eigen::MatrixXd                                                     tipData;
        Eigen::MatrixXd                                                     contrastScatterMatrix;
        Eigen::MatrixXd                                                     varianceCovarianceMatrix;
        Eigen::MatrixXd                                                     psi;
        Eigen::MatrixXd                                                     psiN;
        Eigen::MatrixXd                                                     L;
        Eigen::MatrixXd                                                     C;
        Eigen::MatrixXd                                                     Z;
        Eigen::VectorXd                                                     mu0;
        Eigen::MatrixXd*                                                    originalData;
        ParameterTree*                                                      treeParam;
        Node*                                                               root;
        double                                                              cachedLnL;
        double                                                              cachedLnP;
        bool                                                                branchLengthsInstantiated;
        bool                                                                fixedTree;
        double                                                              dof;
        double                                                              dofN;
        double                                                              lambda;
        double                                                              lambdaN;
        double                                                              log2Pi;
        int                                                                 numberOfInternalNodes;
        int                                                                 numberOfNodes;
        int                                                                 numberOfTips;
        int                                                                 numberOfTraits;
};

#endif /* MultivariateBrownianMotionV2_hpp */
