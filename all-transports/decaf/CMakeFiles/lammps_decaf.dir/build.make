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
include all-transports/decaf/CMakeFiles/lammps_decaf.dir/depend.make

# Include the progress variables for this target.
include all-transports/decaf/CMakeFiles/lammps_decaf.dir/progress.make

# Include the compile flags for this target's objects.
include all-transports/decaf/CMakeFiles/lammps_decaf.dir/flags.make

all-transports/decaf/CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.o: all-transports/decaf/CMakeFiles/lammps_decaf.dir/flags.make
all-transports/decaf/CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.o: all-transports/decaf/lammps_decaf.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object all-transports/decaf/CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.o"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf && /opt/cray/craype/2.5.13/bin/CC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.o -c /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf/lammps_decaf.cpp

all-transports/decaf/CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.i"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf && /opt/cray/craype/2.5.13/bin/CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf/lammps_decaf.cpp > CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.i

all-transports/decaf/CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.s"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf && /opt/cray/craype/2.5.13/bin/CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf/lammps_decaf.cpp -o CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.s

# Object files for target lammps_decaf
lammps_decaf_OBJECTS = \
"CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.o"

# External object files for target lammps_decaf
lammps_decaf_EXTERNAL_OBJECTS =

bin/lammps_decaf: all-transports/decaf/CMakeFiles/lammps_decaf.dir/lammps_decaf.cpp.o
bin/lammps_decaf: all-transports/decaf/CMakeFiles/lammps_decaf.dir/build.make
bin/lammps_decaf: /lustre/atlas/proj-shared/csc143/dhuang/lammps-22Aug18/build/liblammps.a
bin/lammps_decaf: msd-anal/libmsd_analysis.a
bin/lammps_decaf: /usr/lib64/libpng.so
bin/lammps_decaf: /usr/lib64/libjpeg.so
bin/lammps_decaf: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdecaf.a
bin/lammps_decaf: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libmanala.a
bin/lammps_decaf: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_datamodel.a
bin/lammps_decaf: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_transport_mpi.a
bin/lammps_decaf: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbca.a
bin/lammps_decaf: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdca.a
bin/lammps_decaf: all-transports/decaf/CMakeFiles/lammps_decaf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/lammps_decaf"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lammps_decaf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
all-transports/decaf/CMakeFiles/lammps_decaf.dir/build: bin/lammps_decaf

.PHONY : all-transports/decaf/CMakeFiles/lammps_decaf.dir/build

all-transports/decaf/CMakeFiles/lammps_decaf.dir/clean:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf && $(CMAKE_COMMAND) -P CMakeFiles/lammps_decaf.dir/cmake_clean.cmake
.PHONY : all-transports/decaf/CMakeFiles/lammps_decaf.dir/clean

all-transports/decaf/CMakeFiles/lammps_decaf.dir/depend:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/decaf/CMakeFiles/lammps_decaf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : all-transports/decaf/CMakeFiles/lammps_decaf.dir/depend
