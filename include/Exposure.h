#ifndef EXPOSURE_H
#define EXPOSURE_H
#include "RiskMetrics.h"
#include "OptionPricer.h"
#include <vector>

class Exposure {
    public:
    // Define constructor
    Exposure(OptionPricer& option_, double NoOfSim_, double t_);
    // Function
    virtual vector<double> operator()() const;
    virtual double get_t() const;
    virtual void setter_t(double new_t);
    // Define deconstructor
    virtual ~Exposure(){};
    private :
    OptionPricer& option;
    double NoOfSim;
    double t;


    

};


#endif // EXPOSURE_H