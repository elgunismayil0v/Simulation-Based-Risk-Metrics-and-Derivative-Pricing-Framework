cmake_minimum_required(VERSION 3.10)
project(LongstaffSchwartzLSM)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include Eigen3 directory (scoped to this target)
add_executable(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework
    src/main.cpp
    src/Heston_Model_with_Monte_Carlo_simulation.cpp
    src/lmm.cpp
    src/Longstaff_Schwartz_LSM.cpp
)

# Add header directories
target_include_directories(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    /Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/include/eigen-3.4.0/
)

# Add compiler warning flags
target_compile_options(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    -Wall -Wextra
)
