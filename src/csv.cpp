#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Function to write simulation paths to a CSV file (replaces content)
void writeToCSV(const std::vector<std::vector<double> >& data, const std::string& filename) {
    // Open the file in truncation mode (default for ofstream)
    std::ofstream file(filename);
    
    if (file.is_open()) {
        // Write each path (row) to the file
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i < row.size() - 1) file << ",";  // Add comma except after the last value
            }
            file << "\n";  // End the row
        }
        file.close();
        std::cout << "Data written to " << filename << " (file replaced).\n";
    } else {
        std::cerr << "Error: Unable to open file for writing: " << filename << "\n";
    }
}

