#ifndef GAMMA_H
#define GAMMA_H
#include "Greeks.h"
#include "OptionPricer.h"

class Gamma : public Greek {
    public:
    // Define constructor
    Gamma(OptionPricer& option ,double h_); // the perturbation size h
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~Gamma(){};
    private:
    OptionPricer& option;
    double h;

};



#endif // GAMMA_H