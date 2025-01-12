#ifndef LONGSTAFF_SCHWARTZ_H
#define LONGSTAFF_SCHWARTZ_H

#include <vector>
using namespace std;

class Longstaff_Schwartz_LSM {
    public:
    vector<vector<double> > Paths(const vector<vector<double> >& variance_paths,
    int NoOfPaths, int NoOfSteps, double T, double r, double S_0, double kappa, double gamma, double rho, double vbar, double v0); // Generate stock price paths
    double LSM(vector<vector<double> > &paths, double K, double r, double T);  // Longstaff-Schwartz algorithm
    

};

#endif // LONGSTAFF_SCHWARTZ_H