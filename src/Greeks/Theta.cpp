#include <iostream>
#include "OptionPricer.h"
#include "Theta.h"

using namespace std;

Theta::Theta(OptionPricer& option, double h_)
: option(option), h(h_) {}

double Theta::operator()() const {
    double originalT = option.get_T();  // Get the current stock price
    
    // Perturb the stock price upwards and calculate the price
    option.setter_T(originalT + h);   
    double priceUp = option(); 
    
    // Perturb the stock price downwards and calculate the price
    option.setter_T(originalT - h);   
    double priceDown = option(); 
    
    // Restore the original stock price
    option.setter_T(originalT);  
    
    // Calculate delta as the difference in prices divided by the perturbation size
    double theta = (priceUp - priceDown) / (2 * h * 100);
    
    return theta; // Use central difference
}