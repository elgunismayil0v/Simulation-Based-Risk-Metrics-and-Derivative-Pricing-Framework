#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <stdexcept>
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
    double Option_Price();
    double Greek();
    
    void Cos_method();
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

int main() {
    // Instantiate the HestonModel with sample parameters
    HestonModel model;
    model.S_0 = 100.0;    // Initial stock price
    model.v_bar = 0.04;   // Long-term variance
    model.kappa = 2.0;    // Speed of mean reversion
    model.rho = -0.7;     // Correlation coefficient
    model.gamma = 0.5;    // Volatility of variance
    model.r = 0.03;       // Risk-free rate
    model.v0 = 0.04;      // Initial variance

    // Parameters for the characteristic function
    double tau = 1.0;   // Time to maturity
    double u_start = -10.0, u_end = 10.0; // Range for u
    double u_step = 1.0; // Step size for u

    // Loop over u values and compute the characteristic function
    cout << fixed << setprecision(6); // Set precision for better readability
    cout << "u\tReal Part\tImaginary Part" << endl;
    for (double u = u_start; u <= u_end; u += u_step) {
        complex<double> result = model.Heston_ch_function(u, tau, model.gamma, model.kappa, model.rho, model.v0, model.r);
        cout << u << "\t" << real(result) << "\t" << imag(result) << endl;
    }

    return 0;
}





   