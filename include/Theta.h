#ifndef THETA_H
#define THETA_H
#include "Greeks.h"
#include "OptionPricer.h"

class Theta : public Greek {
    public:
    // Define constructor
    Theta(OptionPricer& option ,double h_); // the perturbation size h
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~Theta(){};
    private:
    OptionPricer& option;
    double h;

};



#endif // THETA_H