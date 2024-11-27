#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <numeric>
#include "fftw3.h"
#include "nlopt.h"

using namespace std;
// Define Heston Class
class HestonModel {
    public:
    // Define parametrs
    double S_0; 
    double v_bar;
    double kappa;
    double rho;
    double gamma;
    double r;
    double v0;
    // Define necesseary function
    double Option_Price(double K, double S_0, double r, double tau, int L, int N);
    double Greek();
    
    complex<double> Cos_method(double c, double d, double a, double b, int k);
    complex<double> Heston_ch_function(double u, double tau, 
    double gamma, double kappa, double rho, double v0, double r);
    void Obtimization_funtion();
};

complex<double> HestonModel :: Heston_ch_function(double u, double tau, 
double gamma, double kappa, double rho, double v0, double r){
    complex<double> i(0.0, 1.0);
    // Calculate D
    complex<double> term1 = pow(kappa - i * rho * gamma * u, 2); // (kappa - gamma * rho * i * u)^2
    complex<double> term2 = (pow(u, 2) + i * u) * pow(gamma, 2); // (u^2 + i * u) * gamma^2
    complex<double> D = sqrt(term1 + term2);
    // Calculate G
    complex<double> term11 = (kappa - i * rho* gamma * u - D); //kappa - i * rho * gamma * u - D
    complex<double> term22 =(kappa - i * rho* gamma * u + D);  //kappa - i * rho * gamma * u + D
    complex<double> G = term11 / term22;
    // Calculate A
    complex<double> A_term1 = r * i*u *tau;
    complex<double> A_term2 = (1.0 - exp(-D * tau)) / (1.0 - G * exp(-D * tau)) *
                          (kappa - complex<double>(0.0, 1.0) * rho * gamma * u - D) *
                          (v0 / pow(gamma, 2.0));

    complex<double> A = A_term1 + A_term2;
    // Calculate C
    complex<double> C_term1 = (kappa - i * rho* gamma * u - D) * tau;
    complex<double> C_term2 = ((1.0 - exp(-D * tau) * G) / (1.0 - G)) * complex<double>(2.0);
    complex<double> C_term3 = kappa * v_bar / pow(gamma, 2);
    complex<double> C = C_term3 * (C_term1 - C_term2);
    // Calculate Characteristic Function
    return exp(A + C); 

};

complex<double> HestonModel ::Cos_method(double c, double d, double a, double b, int k){
    // cos_method: Calculates and returns (chi - psi) * 2 / (b - a)
    const double pi = M_PI; // π constant

    // Calculate χ_k(c, d) (ksi)
    double factor = k * pi / (b - a); // Common factor kπ / (b - a)
    double denom = 1 + factor * factor; // Denominator: 1 + (kπ / (b-a))^2

    // Terms for χ_k (ksi)
    double cos_term = cos(factor * (d - a)) * exp(d) - cos(factor * (c - a)) * exp(c);
    double sin_term = sin(factor * (d - a)) * exp(d) - sin(factor * (c - a)) * exp(c);
    double chi = (1.0 / denom) * (cos_term + factor * sin_term);

    // Calculate ψ_k(c, d) (psi)
    double psi;
    if (k == 0) {
        psi = d - c; // Special case when k = 0
    } else {
        double sin_diff = sin(factor * (d - a)) - sin(factor * (c - a));
        psi = sin_diff * (b - a) / (k * pi);
    }

    // Return (chi - psi) * 2 / (b - a)
    return (chi - psi) * 2 / (b - a);
};


const double pi = M_PI; // π constant

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

double HestonModel::Option_Price(double K, double S_0, double r, double tau, int L, int N) {
    // Initialize parameters
    double x = log(S_0 / K);  // Log-moneyness
    double a = -L + sqrt(tau);           // Start of interval
    double b = L + sqrt(tau);            // End of interval
    double c = 0.0;
    double d = b;
    
    // Generate k-values using linspace
    std::vector<int> k_values(N);
    for (int i = 0; i < N; ++i) {
        k_values[i] = i; // k = 0, 1, 2, ..., N-1
    }

    std::complex<double> summation(0.0, 0.0); // To store summation result

    for (int k : k_values) {
        // Calculate u for each k
        double u = k * pi / (b - a);

        // Characteristic function and COS method contribution
        std::complex<double> exp_factor = std::exp(std::complex<double>(0, u * (x - a)));
        std::complex<double> term = Heston_ch_function(u, tau, gamma, kappa, rho, v0, r) *
                                     Cos_method(c, d, a, b, k) * exp_factor;

        // Add to summation
        summation += term;
    }

    // Discount factor
    double discount_factor = K * exp(-r * tau);

    // Calculate final option price
    double result = discount_factor * real(summation);
    return result;
}





int main() {
    // Instantiate the HestonModel
    HestonModel model;
    model.S_0 = 100.0;    // Initial stock price
    model.v_bar = 0.04;   // Long-term variance
    model.kappa = 2.0;    // Speed of mean reversion
    model.rho = -0.7;     // Correlation coefficient
    model.gamma = 0.5;    // Volatility of variance
    model.r = 0.03;       // Risk-free rate
    model.v0 = 0.04;      // Initial variance

    // Parameters for option pricing
    double K = 100;
    double tau = 0.5; // Time to maturity
    int L = 10;      // Interval parameter
    int N = 128;       // Number of Fourier terms
       // Placeholder for u

    // Call Option_Price
    double price = model.Option_Price(K, model.S_0, model.r, tau, L, N);

    // Output the result
    std::cout << "Option Price: " << price << std::endl;

    return 0;
}

