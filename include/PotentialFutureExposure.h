#ifndef POTENTIAL_FUTURE_EXPOSURE_H
#define POTENTIAL_FUTURE_EXPOSURE_H
#include "RiskMetrics.h"
#include "Exposure.h"

class PotentialFutureExposure : public RiskMetrics {
    public:
    // Define constructor
    PotentialFutureExposure(Exposure& exposure_, double quantile_);
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~PotentialFutureExposure(){};
    private:
    Exposure& exposure;
    double quantile;

};




#endif // POTENTIAL_FUTURE_EXPOSURE_H