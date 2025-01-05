cmake_minimum_required(VERSION 3.10)
project(LongstaffSchwartzLSM)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20) # Use C++23 or adjust to C++20 if preferred
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable with source files
add_executable(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework
    src/main.cpp
    src/Heston_Model_with_Monte_Carlo_simulation.cpp
    src/lmm.cpp
    src/Longstaff_Schwartz_LSM.cpp
    src/csv.cpp
)

# Specify the path to the include directory
target_include_directories(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    ${CMAKE_SOURCE_DIR}/include/eigen-3.4.0 # Relative path for Eigen3
    ${CMAKE_SOURCE_DIR}/include/boost_1_87_0 # Relative path for Boost
    ${CMAKE_SOURCE_DIR}/include            # Include custom headers in the 'include' folder
)

# Add include directory for libomp (headers)
include_directories(/opt/homebrew/opt/libomp/include)

# Add library directory for libomp (library files)
link_directories(/opt/homebrew/opt/libomp/lib)





# Include the directory where Python.h and numpy is located
include_directories(/opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/include/python3.13)


include_directories(/opt/homebrew/lib/python3.13/site-packages/numpy/_core/include)

# Link the Python shared library
target_link_libraries(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE /opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/lib/libpython3.13.dylib)




# Add pybind11 include directory from Homebrew

# Pybind11
find_package(pybind11 CONFIG REQUIRED)
target_link_libraries(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE pybind11::module)
include_directories(/opt/homebrew/opt/pybind11/include/pybind11)

# Find Pybind11
find_package(pybind11 REQUIRED)

# Link Pybind11 with your executable
target_link_libraries(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE pybind11::embed)




# Add compiler warning flags for better diagnostics
target_compile_options(Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework PRIVATE
    -Wall -Wextra -pedantic -Wno-unused-parameter 
)





