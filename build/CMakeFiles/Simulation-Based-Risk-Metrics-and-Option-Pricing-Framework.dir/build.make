# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/build

# Include any dependencies generated for this target.
include CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/flags.make

CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/codegen:
.PHONY : CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/codegen

CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o: CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/flags.make
CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o: /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/src/main.cpp
CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o: CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o -MF CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o.d -o CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o -c /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/src/main.cpp

CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/src/main.cpp > CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.i

CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/src/main.cpp -o CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.s

# Object files for target Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework
Simulation__Based__Risk__Metrics__and__Option__Pricing__Framework_OBJECTS = \
"CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o"

# External object files for target Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework
Simulation__Based__Risk__Metrics__and__Option__Pricing__Framework_EXTERNAL_OBJECTS =

Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework: CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/src/main.cpp.o
Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework: CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/build.make
Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework: libSimulationLib.a
Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework: CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/build: Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework
.PHONY : CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/build

CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/clean

CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/depend:
	cd /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/build /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/build /Users/elgun/desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/build/CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/Simulation-Based-Risk-Metrics-and-Option-Pricing-Framework.dir/depend

