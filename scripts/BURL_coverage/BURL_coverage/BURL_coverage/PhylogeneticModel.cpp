#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "PhylogeneticModel.hpp"


PhylogeneticModel::PhylogeneticModel(void) : updatedParameter(nullptr) {
}
 
PhylogeneticModel::~PhylogeneticModel(void) {

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
