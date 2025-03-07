#include <pybind11/pybind11.h>
#include "EuropeanTypeVanillaOptionPricer.h"

namespace py = pybind11;

PYBIND11_MODULE(european_plain_vanilla, m) {
    py::class_<EuropeanPlainVanilla>(m, "EuropeanPlainVanilla")
        .def(py::init<int, double, double, double, EuropeanPlainVanilla::OptionType, double, double, int, double, double, double, double>())
        .def("__call__", &EuropeanPlainVanilla::operator())
        .def("setter_S0", &EuropeanPlainVanilla::setter_S0)
        .def("setter_T", &EuropeanPlainVanilla::setter_T);
    
    py::enum_<EuropeanPlainVanilla::OptionType>(m, "OptionType")
        .value("Put", EuropeanPlainVanilla::OptionType::Put)
        .value("Call", EuropeanPlainVanilla::OptionType::Call)
        .export_values();
}


