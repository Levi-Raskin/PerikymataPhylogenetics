#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

ParameterTree::ParameterTree(double prob, std::vector<std::string> taxonNames, double lam) :  Parameter(prob, "Tree"), lambda(lam), numRejections(0), numAcceptances(0), useCachedLnP(false){
    trees[0] = new Tree(taxonNames, lambda);
    trees[1] = new Tree(*trees[0]);
}

double ParameterTree::lnProbability(void) {
    if(useCachedLnP == false){
        Tree* t = trees[0];
        cachedLnP = 0.0;
        for (Node* p : t->getDownPassSequence())
            {
            if (p != t->getRoot())
                {
                double v = t->getBranchLength(p, p->getAncestor());
                cachedLnP += Probability::Exponential::lnPdf(lambda, v);
                }
            }
        useCachedLnP = true;
    }
    return cachedLnP;
}

void ParameterTree::print(void) {

    trees[0]->print();
}

void ParameterTree::setTree(Tree* t){
    delete trees[0];
    delete trees[1];
    trees[0] = new Tree(*t);
    trees[1] = trees[0];
}

double ParameterTree::update(void) {
    Msg::error("Tried updating parameter tree, but this program is not set up to infer trees");
    double lnP = 0.0;
    return lnP;
}

void ParameterTree::updateForAcceptance(void) {
    numAcceptances++;
    *(trees[1]) = *(trees[0]);
}

void ParameterTree::updateForRejection(void) {
    numRejections++;
    *(trees[0]) = *(trees[1]);
}
