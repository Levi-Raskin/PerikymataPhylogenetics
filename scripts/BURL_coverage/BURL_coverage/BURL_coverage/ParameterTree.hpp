#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <vector>

#include "Parameter.hpp"
#include "Tree.hpp"

class Node;

class ParameterTree : public Parameter {

    public:
                                    ParameterTree(void) = delete;
                                    ParameterTree(double prob, std::vector<std::string> taxonNames, double lam);
        void                        forceBinary(void) { trees[0]->forceBinary(); trees[1]->forceBinary(); }
        double                      getAcceptanceRatio(void) { return ((double) numAcceptances) /( (double)numAcceptances + (double)numRejections ) ;}
        bool                        getAdaptiveProposalActive(void) { return false; }
        Tree*                       getTree(void) { return trees[0]; }
        double                      lnProbability(void);
        void                        print(void);
        void                        setTree(Tree* t);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        
    private:
        //ordered by menmory footprint
        std::vector<std::string>    outgroup;
        Tree*                       trees[2];
        double                      lambda;
        double                      cachedLnP;
        bool                        useCachedLnP;
        int                         numAcceptances;
        int                         numRejections;
};

#endif
