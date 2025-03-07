cmake_minimum_required(VERSION 3.10)
project(LongstaffSchwartzLSM)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find and include pybind11 (for Python bindings)
# find_package(pybind11 REQUIRED)

# Include Eigen, Boost, and general include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include/eigen-3.4.0
    ${CMAKE_SOURCE_DIR}/include/boost_1_87_0
    ${CMAKE_SOURCE_DIR}/include
    /opt/homebrew/opt/libomp/include  # OpenMP
)

# Find Boost
find_package(Boost REQUIRED)

# Manually set the Boost include directory
include_directories(/opt/homebrew/opt/boost/include)

# Alternative:
include_directories(${Boost_INCLUDE_DIRS})

# Add source files for main executable
set(SOURCE_FILES
    src/main.cpp
    src/PathGenerator/PathGenerator.cpp
    src/Vanilla/EuropeanTypeVanillaOptionPricer.cpp
    src/Barrier/EuropeanTypeBarrierDownInOptionPricer.cpp
    src/Barrier/EuropeanTypeBarrierUpInOptionPricer.cpp
    src/Barrier/EuropeanTypeBarrierUpOutOptionPricer.cpp
    src/Barrier/EuropeanTypeBarrierDownOutOptionPricer.cpp
    src/LookBack/EuropeanTypeLookBackOptionPricer.cpp
    src/Asian/EuropeanTypeAsianOptionPricer.cpp
    src/Greeks/Delta.cpp
    src/Greeks/Rho.cpp
    src/Greeks/Gamma.cpp
    src/Greeks/Theta.cpp
    src/Greeks/Vega.cpp
    src/RiskMetrics/Exposure.cpp
    src/RiskMetrics/ExceptedExposure.cpp
    src/RiskMetrics/PotentialFutureExposure.cpp
    src/RiskMetrics/Cva.cpp
    
)

# Create a static/shared library for modularity
add_library(SimulationLib STATIC ${SOURCE_FILES})

# Create the main executable
add_executable(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework src/main.cpp)
target_link_libraries(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE SimulationLib)

# Add the Python module for pybind11 bindings
#pybind11_add_module(european_plain_vanilla src/Bindings/european_plain_vanilla_bindings.cpp)
#target_link_libraries(european_plain_vanilla PRIVATE SimulationLib)

# Link OpenMP (if applicable)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(european_plain_vanilla PRIVATE OpenMP::OpenMP_CXX)
endif()

# Compiler warning flags
target_compile_options(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    -Wall -Wextra -pedantic -Wno-unused-parameter
)
