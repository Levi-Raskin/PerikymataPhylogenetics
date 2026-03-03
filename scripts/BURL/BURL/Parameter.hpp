#ifndef Parameter_hpp
#define Parameter_hpp

#include <string>
class PhylogeneticModel;
class MultivariateNormalPhylogeneticModel;


class Parameter {

    public:
                            Parameter(void) = delete;
                            Parameter(double prob, PhylogeneticModel* m, std::string n);
                            Parameter(double prob, std::string n);
        virtual double      getAcceptanceRatio(void) = 0;
        virtual bool        getAdaptiveProposalActive(void) = 0;
        std::string         getName(void) { return parmName; }
        bool                getParmPrintConsole(void) { return parmPrintsToConsole; }
        double              getProposalProbability(void) { return proposalProbability; }
        virtual double      lnProbability(void) = 0;
        virtual void        print(void) = 0;
        void                setParmPrintConsole(bool x) { parmPrintsToConsole = x; }
        void                setProposalProbability(double x) { proposalProbability = x; }
        virtual double      update(void) = 0;
        virtual void        updateForAcceptance(void) = 0;
        virtual void        updateForRejection(void) = 0;
        
    protected:
        std::string         parmName;
        PhylogeneticModel*  model;
        double              proposalProbability;
        bool                adaptiveProposalActive;
        bool                parmPrintsToConsole;
};

#endif
