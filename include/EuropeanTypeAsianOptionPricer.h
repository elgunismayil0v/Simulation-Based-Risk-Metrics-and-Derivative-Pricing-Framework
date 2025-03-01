#ifndef ASIAN_OPTION_H
#define ASIAN_OPTION_H

#include "EuropeanOption.h"

class AsianOption : public EuropeanOption {
     public:
    enum OptionType {
        Call,
        Put
    };
    enum StrikeType  {
        Fixed,
        Floating
    };
    // Define constructor
    AsianOption(int NoOfPaths_, double strike_, double T_, double r_, OptionType OptionType_, StrikeType StrikeType_, double S0_, double rho_, int NoOfSteps_,
               double kappa_, double gamma_, double vbar_, double v0_);
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~AsianOption(){};

    private:
    int NoOfPaths;
    double strike;
    double T;
    double r;
    OptionType Option_Type;
    StrikeType Strike_Type;
    double S0;
    int NoOfSteps;
    double kappa;
    double gamma;
    double vbar;
    double rho;
    double v0;

};



#endif // ASIAN_OPTION_H
