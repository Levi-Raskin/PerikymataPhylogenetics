#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <vector>

#include "Parameter.hpp"
#include "Tree.hpp"

class Node;

class ParameterTree : public Parameter {

    public:
                                    ParameterTree(void) = delete;
                                    ParameterTree(double prob, PhylogeneticModel* m, std::vector<std::string> taxonNames, std::vector<std::string> og, double lam);
                                    ParameterTree(double prob, PhylogeneticModel* m, std::vector<std::string> taxonNames, double lam);
                                    ParameterTree(double prob, PhylogeneticModel* m, Tree* t, double lam);
        double                      getAcceptanceRatio(void) { return ((double) numAcceptances) /( (double)numAcceptances + (double)numRejections ) ;}
        bool                        getAdaptiveProposalActive(void) { return false; }
        Tree*                       getTree(void) { return trees[0]; }
        double                      lnProbability(void);
        void                        print(void);
        void                        setTree(Tree* t) { trees[0] = t; trees[1] = t; }
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        
    private:
        //Functions
        double                      updateBranchLength(void);
        double                      updateSPR(void);
        //Objects
        std::vector<std::string>    outgroup;
        double                      lambda;
        int                         numAcceptances;
        int                         numRejections;
        Tree*                       trees[2];
};

#endif
