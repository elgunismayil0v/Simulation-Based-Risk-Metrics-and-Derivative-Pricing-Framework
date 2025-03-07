#ifndef CVA_H
#define CVA_H
#include "RiskMetrics.h"
#include "ExceptedExposure.h"
#include "Exposure.h"
#include "OptionPricer.h"

class Cva : public RiskMetrics {
    public:
    // Define constructor
    Cva(Exposure& exposure_, OptionPricer& option_, double hazard_rate_, double recovery_rate_);
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~Cva(){};
    private:
    Exposure& exposure;
    OptionPricer& option;
    double hazard_rate;
    double recovery_rate;



};


#endif // CVA_H