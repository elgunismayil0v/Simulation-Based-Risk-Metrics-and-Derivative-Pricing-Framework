#ifndef DELTA_H
#define DELTA_H
#include "Greeks.h"
#include "OptionPricer.h"

class Delta : public Greek {
    public:
    // Define constructor
    Delta(OptionPricer& option ,double h_); // the perturbation size h
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~Delta(){};
    private:
    OptionPricer& option;
    double h;

};



#endif // DELTA_H