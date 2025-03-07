#ifndef EUROPEAN_PLAIN_VANILLA_H
#define EUROPEAN_PLAIN_VANILLA_H

#include "OptionPricer.h"
#include "EuropeanOption.h"
#include <vector>
using namespace std;



class EuropeanPlainVanilla : public EuropeanOption {
    public:
    enum OptionType {
        Put,
        Call
    };
    // Define constructor
    EuropeanPlainVanilla(int NoOfPaths_ ,double strike_, double T_, double r_, OptionType OptionType_, double S0_, double rho_, int NoOfSteps_,
            double kappa_, double gamma_, double vbar_, double v0_);
    // Function
    virtual double operator()() const override;

    // Virtual getter functions
    virtual double get_T() const override;
    virtual double get_r() const override;
    virtual double get_S0() const override;
    virtual double get_v0() const override;

    // Virtual setter functions
    virtual void setter_T(double new_T_) override;
    virtual void setter_r(double new_r_) override;
    virtual void setter_S0(double new_S0_) override;
    virtual void setter_v0(double new_v0_) override;


    // Define deconstructor
    virtual ~EuropeanPlainVanilla(){};
    private:
    int NoOfPaths;
    double strike;
    double T;
    double r;
    OptionType Option_Type;
    double S0;
    int NoOfSteps;
    double kappa;
    double gamma;
    double vbar;
    double rho;
    double v0;
};

#endif //EUROPEAN_PLAIN_VANILLA_H