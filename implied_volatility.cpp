#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

// Our project has 3 parts. First, we calculate option value use black-scholes.
// Second, derive Vega function (sensitivity of price with respect to volatility).
// Third, Implied volatility using Newton-Raphson method.

//Define functions

struct implied_volatility
{
    double S_0;
    int K;
    double T;
    double interest;
};


// Define cumulative normal dist
double norm_cdf(double x) {
    return 0.5 * erfc(-x * M_SQRT1_2);
}




double black_scholes(const implied_volatility& params, double sigma ,const string& type_of_option) {

    double d1 = (log(params.S_0 / params.K) + (0.5 * sigma * sigma + params.interest) * params.T) / (sigma * sqrt(params.T));

    double d2 = d1 - sigma * sqrt(params.T) ;

    if (type_of_option == "Call"){
        return params.S_0 * norm_cdf(d1) - params.K * norm_cdf(d2) * exp(-params.interest * params.T) ;
    } else if (type_of_option == "Put")
    {
        return params.K * norm_cdf(-d2) * exp(-params.interest * params.T) - params.S_0 * norm_cdf(-d1);
    }

    return 0;
    
}

double Vega(const implied_volatility& params, double sigma) {
    double d1 = (log(params.S_0 / params.K) + (0.5 * sigma * sigma + params.interest) * params.T) / (sigma * sqrt(params.T));
    return params.S_0 * sqrt(params.T) * exp(-0.5 * d1 * d1) / sqrt(2 * M_PI);

}

double Newton_Raphson(const implied_volatility& params, const string& type_of_option, double market_value, int max_iter = 100, double tol = 1e-6) {
     double sigma = 0.2;  // Initial guess for volatility

    for (int i = 0; i < max_iter; i++) {
        double fair_value = black_scholes(params, sigma, type_of_option);
        double diff = market_value - fair_value;

        if (fabs(diff) < tol) return sigma;  // Converged

        double vega_value = Vega(params, sigma);
        if (vega_value < 1e-16) throw runtime_error("Vega too small, Newton-Raphson method failed.");

        sigma -= diff / vega_value;  // Update sigma
    }
    throw runtime_error("Newton-Raphson did not converge within the maximum number of iterations.");
}





int main() {
    implied_volatility params;
    double market_value;
    string type_of_option;
    cout << "Write type of option \n";
    cin >> type_of_option;
    cout << "Write Stock price \n";
    cin >> params.S_0;
    cout << "Write Strike price \n";
    cin >> params.K;
    cout << "Write maturity of Stock \n";
    cin >> params.T;
    cout << "Write interest rate \n";
    cin >> params.interest;
    cout << "Write market value \n";
    cin >> market_value;

    try {
        double implied_vol = Newton_Raphson(params, type_of_option, market_value);
        cout << "Implied Volatility: " << implied_vol << endl;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
    }

    return 0;
    
}

