#ifndef XVA_METRICS_H
#define XVA_METRICS_H

#include <vector>
using namespace std;

class RiskMetrics {
    public:
    // Define constructor
    RiskMetrics(){};
    // Function
    virtual double operator()() const = 0;
    // Define deconstructor
    virtual ~RiskMetrics(){};

};










#endif
