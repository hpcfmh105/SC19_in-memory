# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

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
CMAKE_COMMAND = /autofs/nccs-svm1_sw/titan/.swci/0-login/opt/spack/20180315/linux-suse_linux11-x86_64/gcc-5.3.0/cmake-3.11.3-f2xxt5vasvlzrefneolaxu4rgko7cjnb/bin/cmake

# The command to remove a file.
RM = /autofs/nccs-svm1_sw/titan/.swci/0-login/opt/spack/20180315/linux-suse_linux11-x86_64/gcc-5.3.0/cmake-3.11.3-f2xxt5vasvlzrefneolaxu4rgko7cjnb/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master

# Include any dependencies generated for this target.
include tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/depend.make

# Include the progress variables for this target.
include tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/progress.make

# Include the compile flags for this target's objects.
include tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/flags.make

tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.o: tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/flags.make
tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.o: tests/test-lammps/test_lammps_mpi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.o"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps && /opt/cray/craype/2.5.13/bin/CC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.o -c /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps/test_lammps_mpi.cpp

tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.i"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps && /opt/cray/craype/2.5.13/bin/CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps/test_lammps_mpi.cpp > CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.i

tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.s"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps && /opt/cray/craype/2.5.13/bin/CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps/test_lammps_mpi.cpp -o CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.s

# Object files for target test_lammps_mpi
test_lammps_mpi_OBJECTS = \
"CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.o"

# External object files for target test_lammps_mpi
test_lammps_mpi_EXTERNAL_OBJECTS =

bin/test_lammps_mpi: tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/test_lammps_mpi.cpp.o
bin/test_lammps_mpi: tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/build.make
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/lammps-22Aug18/build/liblammps.a
bin/test_lammps_mpi: /usr/lib64/libpng.so
bin/test_lammps_mpi: /usr/lib64/libjpeg.so
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdecaf.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libmanala.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_datamodel.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_transport_mpi.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbca.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdca.a
bin/test_lammps_mpi: /usr/lib64/libpng.so
bin/test_lammps_mpi: /usr/lib64/libjpeg.so
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdecaf.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libmanala.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_datamodel.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_transport_mpi.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbca.a
bin/test_lammps_mpi: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdca.a
bin/test_lammps_mpi: tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/test_lammps_mpi"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_lammps_mpi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/build: bin/test_lammps_mpi

.PHONY : tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/build

tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/clean:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps && $(CMAKE_COMMAND) -P CMakeFiles/test_lammps_mpi.dir/cmake_clean.cmake
.PHONY : tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/clean

tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/depend:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/depend

