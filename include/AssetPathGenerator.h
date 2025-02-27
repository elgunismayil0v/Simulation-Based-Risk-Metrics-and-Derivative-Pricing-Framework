#ifndef ASSET_GENERATOR_H
#define ASSET_GENERATOR_H

#include <vector>
using namespace std;

class GenerateAssetPaths {
public:
    // Inline default constructor
    GenerateAssetPaths() {}
    // Pure virtual function
    virtual vector<double> operator()(const vector<double>& VarianceProcess) const = 0;
    // Virtual destructor
    virtual ~GenerateAssetPaths() {}
};

class AssetPaths : public GenerateAssetPaths {
public:
    // Define constructor.
    AssetPaths(double S0_, double rho_, double r_, int NoOfSteps_, double T_,
               double kappa_, double gamma_, double vbar_);
    // Compute asset paths based on variance paths.
    virtual vector<double> operator()(const vector<double>& VarianceProcess) const override;
    // Virtual destructor
    virtual ~AssetPaths() {}
private:
    double S0;
    int NoOfSteps;
    double T;
    double kappa;
    double gamma;
    double vbar;
    double rho;
    double r;
};

#endif
