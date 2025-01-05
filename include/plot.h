#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include <string>
#include <vector>


/**
 * @brief Plots a series of stock price paths using matplotlib-cpp.
 * 
 * @param paths A vector of vectors where each inner vector represents a price path.
 * @param title The title of the plot.
 */
void plot_price_paths(const std::vector<std::vector<double> >& paths, const std::string& plot_title);

#endif // PLOT_H
