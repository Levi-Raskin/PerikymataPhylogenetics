#ifndef ParameterDouble_hpp
#define ParameterDouble_hpp

#include "Parameter.hpp"
#include <deque>
#include <string>
#include <vector>

class PhylogeneticModel;

class ParameterDouble : public Parameter {

    public:
                                    ParameterDouble(void) = delete;
                                    ParameterDouble(double prob, PhylogeneticModel* p, std::string n, double min); //min bound
                                    ParameterDouble(double prob, PhylogeneticModel* p, std::string n); //no bounds
        void                        decreaseWindowSize(void) { windowSize *= 0.8; }
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances)/((double)(numAcceptances+numRejections));}
        bool                            getAdaptiveProposalActive(void) { return adaptiveProposalActive; }
        double                      getValue(void) { return value[0]; } // 0 is the one we update, 1 is the one we don't (last currently accepted value
        void                        increaseWindowSize(void) { windowSize *= 1.2; }
        double                      lnProbability(void);
        void                        multiply(double x) { value[1] *= x; } // used in the specific case where we don't have a true parameter, but need to update and be consistent across multiple classes such as mvBM/NIW degree of freedom adaptive tuning
        void                        print(void);
        void                        resetAdaptiveSampling(void) { numRejections = 0; numAcceptances = 0;}
        void                        setNumAdaptive(int x) { numAdaptive = x;}
        void                        setValue(double x) { value[0] = x; value[1] = x;}
        double                      update(void);
        double                      updateAdaptive(int numGen, double targetR);
        double                      updateExp(int numGen, double targetR);
        double                      updateNormal(int numGen, double targetR);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        updateForRejectionNotDynamic(void);
    private:
        int                         numAcceptances;
        int                         numAdaptive;
        int                         numRejections;
        double                      lowerBound;
        double                      upperBound;
        std::vector<double>         value;
        double                      windowSize;
        std::deque<bool>            recentAcceptRej;
};


#endif /* ParameterDouble_hpp */
