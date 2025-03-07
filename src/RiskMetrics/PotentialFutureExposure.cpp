#include <vector>
#include <cmath>  // for exp, sqrt, log
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include "Exposure.h"
#include "PotentialFutureExposure.h"

PotentialFutureExposure::PotentialFutureExposure(Exposure& exposure_, double quantile_)
: exposure(exposure_) {
    quantile = quantile_;
}

#include <algorithm> // For std::sort
#include <vector>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

double PotentialFutureExposure::operator()() const {
    vector<double> sortedExposure = exposure(); 
    sort(sortedExposure.begin(), sortedExposure.end());


    // Corrected index calculation
    double exactIndex = (sortedExposure.size() - 1) * quantile; 
    int Index = static_cast<int>(ceil(exactIndex));


    return sortedExposure[Index];
}



