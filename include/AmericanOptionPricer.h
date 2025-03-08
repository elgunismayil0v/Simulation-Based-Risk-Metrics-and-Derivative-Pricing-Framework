#ifndef AMERICAN_OPTION_H
#define AMERICAN_OPTION_H
#include "OptionPricer.h"

#include <vector>
using namespace std;

class LongstaffSchwartzLSM : public OptionPricer {
public:
     enum OptionType {
        Put,
        Call
    };
    // Define constructor
    LongstaffSchwartzLSM(int NoOfPaths_ ,double strike_, double T_, double r_, OptionType OptionType_, double S0_, double rho_, int NoOfSteps_,
            double kappa_, double gamma_, double vbar_, double v0_);
    // Function
    virtual double operator()() const override;
    vector<vector<double>> matrix() const;
    vector<int> ITM_Paths(const vector<vector<double>>& paths_matrix, int t) const;
    vector<double> payoff(const vector<vector<double>>& paths_matrix) const;

    // Implement setter and getter functions for S0, T, r, v0
    void setter_S0(double new_S0_) override {
        S0 = new_S0_;
    }

    double get_S0() const override {
        return S0;
    }

    void setter_T(double new_T) override {
        T = new_T;
    }

    double get_T() const override {
        return T;
    }

    void setter_r(double new_r) override {
        r = new_r;
    }

    double get_r() const override {
        return r;
    }

    void setter_v0(double new_v0) override {
        v0 = new_v0;
    }

    double get_v0() const override {
        return v0;
    }



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

#endif // AMERICAN_OPTION_H