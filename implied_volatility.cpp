#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

// Our project has 3 parts. First, we calculate option value use black-scholes.
// Second, derive Vega function (sensitivity of price with respect to volatility).
// Third, Implied volatility using Newton-Raphson method.

//Define functions
double black_scholes(double S_0, int K, double d1, double d2, string type_of_option, double sigma,double interest = 0.05,
double value_of_option);

double Vega();

double Newton_Raphson();


double black_scholes(double S_0, int K, double d1, double d2, 
string type_of_option, double sigma,double interest = 0.05, double value_of_option) {
    cout << "Write type of option \n";
    cin >> type_of_option;
    cout << "Write Stock price";
    cin >> S_0;
    cout << "Write Strike price";
    cin >> K;
    cout << "Write volatility of Stock";
    cin >> sigma;

    if (type_of_option == "Call");



}

