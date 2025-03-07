#ifndef OPTION_PRICER_H
#define OPTION_PRICER_H
#include <vector>
using namespace std;



class OptionPricer {
    public:
    // Define constructor
    OptionPricer(){};
    // Functions
    virtual double operator()() const = 0;
    virtual void setter_S0(double new_S0_) = 0;
    virtual double get_S0() const = 0;
    virtual void setter_T(double new_T) = 0;
    virtual double get_T() const = 0;
    virtual void setter_r(double new_r) = 0;
    virtual double get_r() const = 0;
    virtual void setter_v0(double new_v0) = 0;
    virtual double get_v0() const = 0;

    // Define deconstructor
    virtual ~OptionPricer(){};

};

#endif // OPTION_PRICER_H
