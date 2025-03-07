#include <vector>
#include <cmath>  // for exp, sqrt, log
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include "Exposure.h"
#include "OptionPricer.h"
#include "ExceptedExposure.h"

ExceptedExposure :: ExceptedExposure(Exposure& exposure_) : exposure(exposure_) {}

double ExceptedExposure::operator()() const {
    vector<double> values = exposure();  // Retrieve exposures once
    return accumulate(values.begin(), values.end(), 0.0) / values.size();
}
