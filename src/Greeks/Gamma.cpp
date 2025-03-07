#include <iostream>
#include "OptionPricer.h"
#include "Gamma.h"

using namespace std;

Gamma::Gamma(OptionPricer& option, double h_)
: option(option), h(h_) {}

double Gamma::operator()() const {
    double originalS0 = option.get_S0();  // Get the current stock price
    
    // Perturb the stock price upwards and calculate the price
    option.setter_S0(originalS0 + h);   
    double priceUp = option(); 
    
    // Perturb the stock price downwards and calculate the price
    option.setter_S0(originalS0 - h);   
    double priceDown = option(); 
    

    
    // Calculate delta as the difference in prices divided by the perturbation size
    double gamma = (priceUp - 2 * option() + priceDown) / (h * h * 10000);

    // Restore the original stock price
    option.setter_S0(originalS0);  
    
    return gamma;  // Use central difference
}