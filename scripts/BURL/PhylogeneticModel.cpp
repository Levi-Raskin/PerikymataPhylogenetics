#include <iomanip>
#include <iostream>
#include "CharacterMatrix.hpp"
#include "ConditionalLikelihoods.hpp"
#include "DoubleMatrix.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterEigenMatrixD.hpp"
#include "ParameterTree.hpp"
#include "ParameterVarianceCovarianceMatrix.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "StandardMatrix.hpp"
#include "StandardConditionalLikelihoods.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"



PhylogeneticModel::PhylogeneticModel(void) : updatedParameter(nullptr) {
}
 
PhylogeneticModel::~PhylogeneticModel(void) {

}

std::vector<std::string> PhylogeneticModel::getHeaderString(void) {

        std::vector<std::string> headers;
    for (Parameter* p : parameters){
            ParameterVarianceCovarianceMatrix* mat = dynamic_cast<ParameterVarianceCovarianceMatrix*>(p);
            if(mat != nullptr){
                Eigen::MatrixXd scratch = (mat->getValue());
                for(int i = 0; i < scratch.rows(); i++)
                    for(int j = 0; j < scratch.cols(); j++)
                        headers.push_back( p->getName() + std::to_string(i) + "," + std::to_string(j));
            }else{
                headers.push_back( p->getName() );
            }
        }
//    headers.push_back("adaptiveSamplingDone(0false1true)");
    return headers;
}

Tree* PhylogeneticModel::getTree(void) {
    for (Parameter* p : parameters)
        {
        ParameterTree* pt = dynamic_cast<ParameterTree*>(p);
        if (pt != nullptr)
            return pt->getTree();
        }
    return nullptr;
}
