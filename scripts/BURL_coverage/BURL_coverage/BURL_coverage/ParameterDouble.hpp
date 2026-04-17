#ifndef ParameterDouble_hpp
#define ParameterDouble_hpp

#include "Parameter.hpp"
#include <deque>
#include <string>
#include <vector>

class ParameterDouble : public Parameter {

    public:
                                    ParameterDouble(void) = delete;
                                    ParameterDouble(double prob, std::string n);
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances)/((double)(numAcceptances+numRejections));}
        bool                            getAdaptiveProposalActive(void) { return adaptiveProposalActive; }
        double                      getValue(void) { return value[0]; } // 0 is the one we update, 1 is the one we don't (last currently accepted value
        double                      lnProbability(void);
        void                        print(void);
        void                        setValue(double x) { value[0] = x;}
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        updateForRejectionNotDynamic(void);
    private:
        double                      updateAdaptive(int numGen, double targetR);
        std::vector<double>         value;
        std::deque<bool>            recentAcceptRej;
        double                      lowerBound;
        double                      upperBound;
        double                      windowSize;
        int                         numAcceptances;
        int                         numAdaptive;
        int                         numRejections;
};


#endif /* ParameterDouble_hpp */
