#ifndef VARIANCE_GENERATOR_H
#define VARIANCE_GENERATOR_H

#include <vector>
using namespace std;


//Base class for variance process generators.
class GenerateVariancePaths {
public:
    GenerateVariancePaths(){}
    virtual vector<double> operator()() const = 0;
    virtual ~GenerateVariancePaths(){}
};

//Heston Variance Process.
class VarianceProcess : public GenerateVariancePaths {
public:
    VarianceProcess(int NoOfSteps_, double T_,
                    double kappa_, double gamma_, double vbar_, double v0_);
    
    virtual vector<double> operator()() const override;
    virtual ~VarianceProcess(){};
    
private:
    int NoOfSteps;
    double T;
    double kappa;
    double gamma;
    double vbar;
    double v0;
};

#endif // VARIANCE_GENERATOR_H
