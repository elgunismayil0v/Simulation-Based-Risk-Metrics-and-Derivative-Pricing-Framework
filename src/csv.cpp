#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Function to append data (simulation paths) to a CSV file
void appendToCSV(const std::vector<std::vector<double> >& data, const std::string& filename) {
    // Open the file in append mode
    std::ofstream file(filename, std::ios::app);
    
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
        std::cout << "Data appended to " << filename << "\n";
    } else {
        std::cerr << "Error: Unable to open file for writing: " << filename << "\n";
    }
}
