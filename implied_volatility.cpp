#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

// Our project has 3 parts. First, we calculate option value use black-scholes.
// Second, derive Vega function (sensitivity of price with respect to volatility).
// Third, Implied volatility using Newton-Raphson method.

//Define functions
double black_scholes(double S_0, int K, double d1, double d2, string type_of_option, double sigma,double interest ,
double value_of_option, double T);

double Vega();

double Newton_Raphson();

// Define cumulative normal dist
double norm_cdf(double x) {
    return 0.5 * erfc(-x * M_SQRT1_2);
}





double black_scholes(double S_0, int K,  
string type_of_option, double sigma, double interest, double value_of_option, double T) {
    

    double d1 = (log(S_0/K) + (0.5 * sigma * sigma + interest) * T) / (sigma * sqrt(T)) ;
    double d2 = d1 - sigma * sqrt(T) ;


    if (type_of_option == "Call"){
        return  S_0 * norm_cdf(d1) - K * norm_cdf(d2) * exp(-interest * T) ;
    } else if (type_of_option == "Put")
    {
        return  K * norm_cdf(-d2) * exp(-interest * T) - S_0 * norm_cdf(-d1);
    }

   

    return 0;
    
}

int main() {
    double sigma , S_0, K, T, interest;
    string type_of_option;
    cout << "Write type of option \n";
    cin >> type_of_option;
    cout << "Write Stock price";
    cin >> S_0;
    cout << "Write Strike price";
    cin >> K;
    cout << "Write volatility of Stock";
    cin >> sigma;
    cout << "Write time to maturity of option";
    cin >> T;
    cout << "Write interest rate";
    cin >> interest;

    cin.get() ;
    cin.get() ;
    
}

