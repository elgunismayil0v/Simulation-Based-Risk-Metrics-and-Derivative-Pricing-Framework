#include "EuropeanOption.h"

using namespace std;

void EuropeanOption::setter_S0( double new_S0_) {
    S0 = new_S0_;
}

double EuropeanOption::get_S0() const {
    return S0; 
}

void EuropeanOption::setter_r( double new_r_) {
    r = new_r_;
}

double EuropeanOption::get_r() const {
    return r; 
}

void EuropeanOption::setter_v0( double new_v0_) {
    v0 = new_v0_;
}

double EuropeanOption::get_v0() const {
    return v0; 
}

void EuropeanOption::setter_T( double new_T_) {
    T = new_T_;
}

double EuropeanOption::get_T() const {
    return T; 
}