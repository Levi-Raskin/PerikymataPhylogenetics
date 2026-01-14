#ifndef MultivariateBrownianMotion_hpp
#define MultivariateBrownianMotion_hpp

#include "Eigen/Dense"
#include "PhylogeneticModel.hpp"
#include <map>
#include <string>
#include <vector>

class Tree;
class Node;
class Parameter;
class ParameterDouble;
class ParameterMatrix;
class ParameterVectorDouble;

class MultivariateBrownianMotion : public PhylogeneticModel{
    public:
                                                                            MultivariateBrownianMotion(std::vector<std::string> rn, Eigen::MatrixXd* data);
        std::vector<std::string>                                            getParameterNames(void);
        std::vector<double>                                                 getParameterString(void);
        double                                                              lnLikelihood(void);
        double                                                              lnPriorProbability(void);
        void                                                                print(void);
        void                                                                simulateData(void);
        void                                                                setTree(Tree* t);
        double                                                              update(void);
        void                                                                updateForAcceptance(void);
        void                                                                updateForRejection(void);

    protected:
        void                                                                addData(std::vector<std::string> rn, Eigen::MatrixXd* data);
        std::map<int, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>>          ancDists;
        double                                                              lnLikelihoodF85Optim(void);
        int                                                                 numberOfTraits;
        Eigen::MatrixXd*                                                    observedData;
        Eigen::MatrixXd                                                     datScatter;
        std::vector<std::string>                                            observedDataRownames;
        std::unordered_map<std::string, int>                                tipHashMap;
        std::vector<std::string>                                            tipMeanOrdering;
        ParameterMatrix*                                                    varianceCovarianceMatrix;
};

#endif
