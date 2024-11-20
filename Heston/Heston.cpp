#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <fftw3.h>
#include <nlopt.hpp>

using namespace std;

// Define complex number type for convenience
typedef complex<double> Complex;

// Heston model's characteristic function
Complex hestonCharacteristicFunction(Complex u, double tau, double kappa, double theta, double sigma,
                                      double rho, double v0, double r, double q, double S0) {
    Complex i(0.0, 1.0); // Imaginary unit
    double x0 = log(S0);
    double a = kappa * theta;
    double b = kappa;
    Complex d = sqrt(pow(rho * sigma * u * i - b, 2.0) - sigma * sigma * (u * i + u * u));
    Complex g = (b - rho * sigma * u * i + d) / (b - rho * sigma * u * i - d);

    Complex C = r * u * i * tau + (a / (sigma * sigma)) * ((b - rho * sigma * u * i + d) * tau - 2.0 * log((1.0 - g * exp(d * tau)) / (1.0 - g)));
    Complex D = ((b - rho * sigma * u * i + d) / (sigma * sigma)) * ((1.0 - exp(d * tau)) / (1.0 - g * exp(d * tau)));
    return exp(C + D * v0 + i * u * x0);
}

// FFT-based option pricing using FFTW
void fftOptionPricing(double S0, double tau, double r, double q, double kappa, double theta, double sigma,
                      double rho, double v0, const vector<double>& strikes, vector<double>& modelPrices) {
    const int N = 1024;          // Number of FFT points
    const double alpha = 1.5;   // Damping factor
    const double eta = 0.25;    // Grid spacing in frequency domain
    const double b = M_PI / eta; // Upper bound of integration range

    // Allocate memory for FFTW input and output
    fftw_complex *fftInput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *fftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Create FFTW plan
    fftw_plan plan = fftw_plan_dft_1d(N, fftInput, fftOutput, FFTW_FORWARD, FFTW_ESTIMATE);

    // Populate the input array
    for (int k = 0; k < N; k++) {
        double v = k * eta; // Frequency domain variable
        Complex u = Complex(v, 0) - (alpha + 1.0) * Complex(0, 1); // Ensure v is treated as Complex
        Complex phi = hestonCharacteristicFunction(u, tau, kappa, theta, sigma, rho, v0, r, q, S0);
        Complex numerator = phi * exp(-r * tau) * exp(-alpha * log(S0));
        Complex denominator = Complex(alpha * alpha + alpha - v * v, 2.0 * alpha * v);
        Complex psi = numerator / denominator;

        // Assign real and imaginary parts to fftw_complex
        fftInput[k][0] = real(psi * exp(Complex(0, b * k)) * eta); // Real part
        fftInput[k][1] = imag(psi * exp(Complex(0, b * k)) * eta); // Imaginary part
    }

    // Execute FFT
    fftw_execute(plan);

    // Extract option prices from FFT output
    modelPrices.clear();
    for (const double& strike : strikes) {
        double multiplier = exp(-alpha * log(strike));
        int k = round((log(strike) + b) / eta);
        double price = fftOutput[k][0] * multiplier / M_PI; // Real part of FFT result
        modelPrices.push_back(price);
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(fftInput);
    fftw_free(fftOutput);
}

// CSV Reading Function
void readCSV(const string& fileName, vector<double>& strikes, vector<double>& marketPrices, vector<double>& maturities) {
    ifstream file(fileName);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << fileName << endl;
        return;
    }

    string line;
    getline(file, line); // Skip header
    while (getline(file, line)) {
        stringstream ss(line);
        string strikeStr, priceStr, maturityStr;

        // Parse the line (strike, price, maturity)
        if (getline(ss, strikeStr, ',') && getline(ss, priceStr, ',') && getline(ss, maturityStr, ',')) {
            strikes.push_back(stod(strikeStr));
            marketPrices.push_back(stod(priceStr));
            maturities.push_back(stod(maturityStr));
        }
    }

    file.close();
}

// Objective function for optimization
double objectiveFunction(const vector<double>& params, vector<double> marketPrices, vector<double> strikes,
                         vector<double> maturities, double S0, double r, double q) {
    double kappa = params[0];
    double theta = params[1];
    double sigma = params[2];
    double rho = params[3];
    double v0 = params[4];

    double error = 0.0;

    // Price each option using its maturity
    for (size_t i = 0; i < strikes.size(); i++) {
        vector<double> modelPrices;
        fftOptionPricing(S0, maturities[i], r, q, kappa, theta, sigma, rho, v0, {strikes[i]}, modelPrices);
        error += pow(modelPrices[0] - marketPrices[i], 2);
    }

    return error;
}

// Main function
int main() {
    // Heston model parameters (initial guess)
    double S0 = 4000.0;       // Spot price (example)
    double r = 0.03;          // Risk-free rate
    double q = 0.02;          // Dividend yield
    vector<double> initialParams = {2.0, 0.04, 0.2, -0.7, 0.04};

    // Observed market prices and strikes
    vector<double> strikes, marketPrices, maturities;

    // Read market data from CSV
    string fileName = "data_for_project.csv"; // Replace with your actual file path
    readCSV(fileName, strikes, marketPrices, maturities);

    if (strikes.empty() || marketPrices.empty() || maturities.empty()) {
        cerr << "Error: No data read from CSV file." << endl;
        return -1;
    }

    // Set up optimization
    nlopt::opt opt(nlopt::LN_NELDERMEAD, initialParams.size());
    opt.set_min_objective([](const vector<double>& params, vector<double>& grad, void* data) {
        auto* d = static_cast<tuple<vector<double>*, vector<double>*, vector<double>*>*>(data);
        return objectiveFunction(params, *get<0>(*d), *get<1>(*d), *get<2>(*d), 4000.0, 0.03, 0.02);
    }, new tuple(&marketPrices, &strikes, &maturities));
    opt.set_xtol_rel(1e-6);
    opt.set_lower_bounds({0.01, 0.01, 0.01, -1.0, 0.01});
    opt.set_upper_bounds({10.0, 1.0, 1.0, 1.0, 1.0});

    // Optimize parameters
    double minf;
    nlopt::result result = opt.optimize(initialParams, minf);

    // Output results
    cout << "Optimized Parameters:" << endl;
    cout << "kappa = " << initialParams[0] << endl;
    cout << "theta = " << initialParams[1] << endl;
    cout << "sigma = " << initialParams[2] << endl;
    cout << "rho = " << initialParams[3] << endl;
    cout << "v0 = " << initialParams[4] << endl;

    return 0;
}








