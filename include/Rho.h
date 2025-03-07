#ifndef RHO_H
#define RHO_H
#include "Greeks.h"
#include "OptionPricer.h"

class Rho : public Greek {
    public:
    // Define constructor
    Rho(OptionPricer& option ,double h_); // the perturbation size h
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~Rho(){};
    private:
    OptionPricer& option;
    double h;

};



#endif // DELTA_H