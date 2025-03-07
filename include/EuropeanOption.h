#ifndef EUROPEAN_OPTION_H
#define EUROPEAN_OPTION_H

#include "OptionPricer.h"
#include <vector>
using namespace std;

class EuropeanOption : public OptionPricer {
    public: 
    // Define constructor
    EuropeanOption(){};
    // Function
    virtual double operator()() const = 0;
    
    // Define deconstructor
    virtual ~EuropeanOption(){};

};

#endif //EUROPEAN_OPTION_H