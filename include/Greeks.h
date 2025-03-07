#ifndef GREEKS_H
#define GREEKS_H

class Greek {
    public:
    // Define constructor
    Greek(){};
    // Define function
    virtual double operator()() const = 0;
    // Define deconstructor
    virtual ~Greek(){};
};


#endif // GREEKS_H