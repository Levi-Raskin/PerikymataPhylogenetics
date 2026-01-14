#include "Parameter.hpp"
#include "ParameterTree.hpp"
#include "PhylogeneticModel.hpp"

PhylogeneticModel::PhylogeneticModel(void) {
 
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
