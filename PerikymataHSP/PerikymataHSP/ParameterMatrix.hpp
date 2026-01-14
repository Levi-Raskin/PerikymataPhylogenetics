#ifndef ParameterMatrix_hpp
#define ParameterMatrix_hpp

//base class for all matrix-type parameters, expects eigen matrixxd return type
#include "Eigen/Dense"
#include "Parameter.hpp"

class PhylogeneticModel;

class ParameterMatrix : public Parameter{
    public:
                                ParameterMatrix(double prob, PhylogeneticModel* p, std::string n);
        PhylogeneticModel*      getModel(void) { return model; }
        virtual Eigen::MatrixXd getValue(void) = 0;
        virtual double          lnProbability(void) = 0;
        virtual void            print(void) = 0;
        virtual double          update(void) = 0;
        virtual void            updateForAcceptance(void) = 0;
        virtual void            updateForRejection(void) = 0;
};

#endif
