#include <vector>
#include <iostream>
#include <cmath>  // for exp, sqrt, log
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include "Exposure.h"
#include "Exposure.h"
#include "ExceptedExposure.h"
#include "OptionPricer.h"
#include "Cva.h"

Cva::Cva(Exposure& exposure_, OptionPricer& option_, double hazard_rate_, double recovery_rate_)
: exposure(exposure_), option(option_) {
    hazard_rate = hazard_rate_;
    recovery_rate = recovery_rate_;

}

double Cva::operator()() const {
    double maturity = option.get_T();
    double t_orginal = exposure.get_t();
    double r = option.get_r();
    double dt = 0.25;
    double cva = 0.0;

    for (double t = 0; t < maturity ; t += dt) {
        exposure.setter_t(t);
        vector<double> exposures = exposure();
        double ee = accumulate(exposures.begin(), exposures.end(), 0.0) / exposures.size();
        double S_t = exp(-hazard_rate * t);
        double S_t_dt = exp(-hazard_rate * (t + dt));
        double default_prob = S_t - S_t_dt;  // Default probability in (t, t+dt)
        cva += (1 - recovery_rate) * ee * default_prob * exp(-r * t);
    }

    exposure.setter_t(t_orginal);
    return cva;
}
