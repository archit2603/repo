# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/jinbo/cmake-3.6.0-rc4-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/jinbo/cmake-3.6.0-rc4-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jinbo/CCMpred

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jinbo/CCMpred

# Include any dependencies generated for this target.
include lib/libconjugrad/CMakeFiles/polynomial.dir/depend.make

# Include the progress variables for this target.
include lib/libconjugrad/CMakeFiles/polynomial.dir/progress.make

# Include the compile flags for this target's objects.
include lib/libconjugrad/CMakeFiles/polynomial.dir/flags.make

lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o: lib/libconjugrad/CMakeFiles/polynomial.dir/flags.make
lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o: lib/libconjugrad/samples/polynomial.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jinbo/CCMpred/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o"
	cd /home/jinbo/CCMpred/lib/libconjugrad && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/polynomial.dir/samples/polynomial.c.o   -c /home/jinbo/CCMpred/lib/libconjugrad/samples/polynomial.c

lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/polynomial.dir/samples/polynomial.c.i"
	cd /home/jinbo/CCMpred/lib/libconjugrad && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jinbo/CCMpred/lib/libconjugrad/samples/polynomial.c > CMakeFiles/polynomial.dir/samples/polynomial.c.i

lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/polynomial.dir/samples/polynomial.c.s"
	cd /home/jinbo/CCMpred/lib/libconjugrad && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jinbo/CCMpred/lib/libconjugrad/samples/polynomial.c -o CMakeFiles/polynomial.dir/samples/polynomial.c.s

lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.requires:

.PHONY : lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.requires

lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.provides: lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.requires
	$(MAKE) -f lib/libconjugrad/CMakeFiles/polynomial.dir/build.make lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.provides.build
.PHONY : lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.provides

lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.provides.build: lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o


# Object files for target polynomial
polynomial_OBJECTS = \
"CMakeFiles/polynomial.dir/samples/polynomial.c.o"

# External object files for target polynomial
polynomial_EXTERNAL_OBJECTS =

lib/libconjugrad/bin/polynomial: lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o
lib/libconjugrad/bin/polynomial: lib/libconjugrad/CMakeFiles/polynomial.dir/build.make
lib/libconjugrad/bin/polynomial: lib/libconjugrad/libconjugrad.a
lib/libconjugrad/bin/polynomial: lib/libconjugrad/CMakeFiles/polynomial.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jinbo/CCMpred/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable bin/polynomial"
	cd /home/jinbo/CCMpred/lib/libconjugrad && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/polynomial.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/libconjugrad/CMakeFiles/polynomial.dir/build: lib/libconjugrad/bin/polynomial

.PHONY : lib/libconjugrad/CMakeFiles/polynomial.dir/build

lib/libconjugrad/CMakeFiles/polynomial.dir/requires: lib/libconjugrad/CMakeFiles/polynomial.dir/samples/polynomial.c.o.requires

.PHONY : lib/libconjugrad/CMakeFiles/polynomial.dir/requires

lib/libconjugrad/CMakeFiles/polynomial.dir/clean:
	cd /home/jinbo/CCMpred/lib/libconjugrad && $(CMAKE_COMMAND) -P CMakeFiles/polynomial.dir/cmake_clean.cmake
.PHONY : lib/libconjugrad/CMakeFiles/polynomial.dir/clean

lib/libconjugrad/CMakeFiles/polynomial.dir/depend:
	cd /home/jinbo/CCMpred && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jinbo/CCMpred /home/jinbo/CCMpred/lib/libconjugrad /home/jinbo/CCMpred /home/jinbo/CCMpred/lib/libconjugrad /home/jinbo/CCMpred/lib/libconjugrad/CMakeFiles/polynomial.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/libconjugrad/CMakeFiles/polynomial.dir/depend

