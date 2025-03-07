#include <vector>
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
    for(double t = 0; t < maturity; t += dt) {
        exposure.setter_t(t);
        double ee = accumulate(exposure().begin(), exposure().end(), 0.0) / exposure().size();
        cva += (1 - recovery_rate) * ee * exp(- hazard_rate) * exp(- r * t) * dt;
    }
    exposure.setter_t(t_orginal);
    return cva;
    
}