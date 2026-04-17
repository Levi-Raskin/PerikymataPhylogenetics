#ifndef PhylogeneticModel_hpp
#define PhylogeneticModel_hpp

#include <string>
#include <vector>
class Parameter;
class Tree;

class PhylogeneticModel {

    public:
                                                PhylogeneticModel(void);
                                               ~PhylogeneticModel(void);
        std::vector<std::string>                getHeaderString(void);
        virtual std::vector<std::string>        getParameterNames(void) = 0;
        virtual std::vector<double>             getParameterString(void) = 0;
        Tree*                                   getTree(void);
        const Parameter*                        getUpdatedParameter(void) { return updatedParameter; };
        virtual double                          lnLikelihood(void) = 0;
        virtual double                          lnPriorProbability(void) = 0;
        virtual void                            print(void) = 0;
        virtual double                          update(void) = 0;
        virtual void                            updateForAcceptance(void) = 0;
        virtual void                            updateForRejection(void) = 0;
        std::vector<Parameter*>                 parameters;
        Parameter*                              updatedParameter;
};

#endif
