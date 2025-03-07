#include <vector>
#include <cmath>  // for exp, sqrt, log
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include "Exposure.h"
#include "OptionPricer.h"

Exposure::Exposure(OptionPricer& option_, double NoOfSim_, double t_) : option(option_) {
    NoOfSim = NoOfSim_;
    t = t_;
}

vector<double> Exposure::operator()() const {
    double orginal_T = option.get_T();
    double theta = orginal_T - t;
    double r = option.get_r();
    vector<double>exposure(NoOfSim);
    option.setter_T(theta);
    for (int i = 0; i < NoOfSim; i++) {
        exposure[i] = option() * exp(r * theta);
    }
    option.setter_T(orginal_T);

    return exposure;
    

}

double Exposure::get_t() const {
    return t;

}

void Exposure::setter_t(double new_t)  {
    t = new_t;

}