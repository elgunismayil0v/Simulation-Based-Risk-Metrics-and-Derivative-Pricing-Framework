#ifndef OPTION_PRICER_H
#define OPTION_PRICER_H
#include <vector>
using namespace std;



class OptionPricer {
    public:
    // Define constructor
    OptionPricer(){};
    // Function
    virtual double operator()() const = 0;
    // Define deconstructor
    virtual ~OptionPricer(){};

};

#endif // OPTION_PRICER_H
