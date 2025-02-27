cmake_minimum_required(VERSION 3.10)
project(LongstaffSchwartzLSM)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable with source files
add_executable(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework
    src/main.cpp
    src/Simulation_1.cpp
    src/EuropeanTypeOptionPricer.cpp
    src/EuropeanTypeBarrierDownInOptionPricer.cpp
    src/EuropeanTypeBarrierUpInOptionPricer.cpp
    src/EuropeanTypeBarrierUpOutOptionPricer.cpp
    src/EuropeanTypeBarrierDownOutOptionPricer.cpp
   
)

# Include Eigen and Boost
target_include_directories(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    ${CMAKE_SOURCE_DIR}/include/eigen-3.4.0
    ${CMAKE_SOURCE_DIR}/include/boost_1_87_0
    ${CMAKE_SOURCE_DIR}/include/

)

# Find Python and NumPy
find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
target_include_directories(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
)
target_link_libraries(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE ${Python_LIBRARIES})

# Find and link Pybind11
find_package(pybind11 REQUIRED)
target_link_libraries(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE pybind11::embed)

# Include libomp
include_directories(/opt/homebrew/opt/libomp/include)
link_directories(/opt/homebrew/opt/libomp/lib)

# Add compiler warning flags for better diagnostics
target_compile_options(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    -Wall -Wextra -pedantic -Wno-unused-parameter
)
