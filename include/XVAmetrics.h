#ifndef XVA_METRICS_H
#define XVA_METRICS_H

#include <vector>
using namespace std;

class XVAmetrics {
    public:
    // Define constructor
    XVAmetrics(){};
    // Function
    virtual double operator()(vector<vector<double>> & AssetPaths) const = 0;
    // Define deconstructor
    virtual ~XVAmetrics(){};

};

class ExpectedExposure : public XVAmetrics {
    public:
    // Define constructor
    ExpectedExposure(){};
    // Function
    virtual double operator()(vector<vector<double>> & AssetPaths) const;
    // Define deconstructor
    virtual ~ExpectedExposure(){};

};

class PotentialFutureExposure : public XVAmetrics {

};

class CVA: public XVAmetrics {

};








#endif
