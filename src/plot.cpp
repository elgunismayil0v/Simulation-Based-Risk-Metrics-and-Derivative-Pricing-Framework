#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <iostream>

namespace py = pybind11;

void plot_price_paths(const std::vector<std::vector<double>>& paths, const std::string& plot_title) {
    // Initialize Python interpreter
    py::scoped_interpreter guard{};

    // Import required modules
    py::module_ plt = py::module_::import("matplotlib.pyplot");
    py::module_ np = py::module_::import("numpy");

    // Convert paths to a NumPy array
    py::array_t<double> numpy_array({paths.size(), paths[0].size()});
    for (size_t i = 0; i < paths.size(); ++i) {
        for (size_t j = 0; j < paths[i].size(); ++j) {
            *numpy_array.mutable_data(i, j) = paths[i][j];
        }
    }

    // Plot the data
    plt.attr("figure")();
    for (size_t i = 0; i < paths.size(); ++i) {
        plt.attr("plot")(numpy_array.attr("__getitem__")(i));
    }

    plt.attr("xlabel")("Time Step");
    plt.attr("ylabel")("Stock Price");
    plt.attr("title")(plot_title);
    plt.attr("show")();
}

