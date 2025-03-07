#include <iostream>
#include "OptionPricer.h"
#include "Rho.h"

using namespace std;

Rho::Rho(OptionPricer& option, double h_)
: option(option), h(h_) {}

double Rho::operator()() const {
    double original_r = option.get_r();  // Get the current stock price
    
    // Perturb the stock price upwards and calculate the price
    option.setter_r(original_r + h);   
    double priceUp = option(); 
    
    // Perturb the stock price downwards and calculate the price
    option.setter_r(original_r - h);   
    double priceDown = option(); 
    
    // Restore the original stock price
    option.setter_r(original_r);  
    
    // Calculate delta as the difference in prices divided by the perturbation size
    double rho = (priceUp - priceDown) / (2 * h * 100);
    
    return rho;  // Use central difference
}