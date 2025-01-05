#ifndef CSV_H
#define CSV_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Function to append data (simulation paths) to a CSV file
void appendToCSV(const std::vector<std::vector<double> >& data, const std::string& filename);

#endif // CSV_H
