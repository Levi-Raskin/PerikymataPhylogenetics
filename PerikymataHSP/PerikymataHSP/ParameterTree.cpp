#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

ParameterTree::ParameterTree(double prob, PhylogeneticModel* m, std::vector<std::string> taxonNames, std::vector<std::string> og, double lam) : Parameter(prob, m, "Tree"), lambda(lam), outgroup(og), numRejections(0), numAcceptances(0){

    trees[0] = new Tree(taxonNames, lambda);
    trees[1] = new Tree(*trees[0]);
}

ParameterTree::ParameterTree(double prob, PhylogeneticModel* m, std::vector<std::string> taxonNames, double lam) :  Parameter(prob, m, "Tree"), lambda(lam), numRejections(0), numAcceptances(0){
    trees[0] = new Tree(taxonNames, lambda);
    trees[1] = new Tree(*trees[0]);
}

ParameterTree::ParameterTree(double prob, PhylogeneticModel* m, Tree* t, double lam) :  Parameter(prob, m, "Tree"), lambda(lam), numRejections(0), numAcceptances(0){
    trees[0] = t;
    trees[1] = new Tree(*trees[0]);
}

double ParameterTree::lnProbability(void) {
    Tree* t = trees[0];
    double lnP = 0.0;
    for (Node* p : t->getDownPassSequence()){
        if (p != t->getRoot())
            {
            double v = p->getBranchLength();
            lnP += Probability::Exponential::lnPdf(lambda, v);
            }
        }
    return lnP;
}

void ParameterTree::print(void) {

    trees[0]->print();
}

double ParameterTree::update(void) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double lnP = 0.0;
    if (rng.uniformRv() < 0.5)
        lnP = updateBranchLength();
    else
        lnP = updateSPR();
        
    Msg::error("not working");
    return lnP;
}

double ParameterTree::updateBranchLength(void){
    Msg::error("move not workign yet");
    return -1.0;
}

double ParameterTree::updateSPR(void) {
    Tree* t = trees[0];
    Msg::error("SPR move not workign yet");
    return 0.0;
}

void ParameterTree::updateForAcceptance(void) {
    numAcceptances++;
    *(trees[1]) = *(trees[0]);
}

void ParameterTree::updateForRejection(void) {
    numRejections++;
    *(trees[0]) = *(trees[1]);
}
