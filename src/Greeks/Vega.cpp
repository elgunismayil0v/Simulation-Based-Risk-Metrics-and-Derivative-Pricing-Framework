#include <iostream>
#include "OptionPricer.h"
#include "Vega.h"

using namespace std;

Vega::Vega(OptionPricer& option, double h_)
: option(option), h(h_) {}

double Vega::operator()() const {
    double originalv0 = option.get_v0();  // Get the current stock price
    
    // Perturb the stock price upwards and calculate the price
    option.setter_v0(originalv0 + h);   
    double priceUp = option(); 
    
    // Perturb the stock price downwards and calculate the price
    option.setter_v0(originalv0 - h);   
    double priceDown = option(); 
    
    // Restore the original stock price
    option.setter_v0(originalv0);  
    
    // Calculate delta as the difference in prices divided by the perturbation size
    double vega = (priceUp - priceDown) / (2 * h * 100);
    
    return vega; // Use central difference
}