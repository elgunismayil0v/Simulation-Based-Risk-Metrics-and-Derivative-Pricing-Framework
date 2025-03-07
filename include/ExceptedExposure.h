#ifndef EXCEPTEDEXPOSURE_H
#define EXCEPTEDEXPOSURE_H
#include "Exposure.h"
#include "RiskMetrics.h"

class ExceptedExposure : public RiskMetrics {
    public:
    // Define constructor
    ExceptedExposure(Exposure& exposure_);
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~ExceptedExposure(){};
    private:
    Exposure& exposure;

};



#endif //EXCEPTEDEXPOSURE_H