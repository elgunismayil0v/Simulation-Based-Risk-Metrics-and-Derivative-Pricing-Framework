#ifndef VEGA_H
#define VEGA_H
#include "Greeks.h"
#include "OptionPricer.h"

class Vega : public Greek {
    public:
    // Define constructor
    Vega(OptionPricer& option ,double h_); // the perturbation size h
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~Vega(){};
    private:
    OptionPricer& option;
    double h;

};



#endif // VEGA_H