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
include all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/depend.make

# Include the progress variables for this target.
include all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/progress.make

# Include the compile flags for this target's objects.
include all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/flags.make

all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.o: all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/flags.make
all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.o: all-transports/adios-all/heat3d-adios/heat3d_mpiio_con.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.o"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.o   -c /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios/heat3d_mpiio_con.c

all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.i"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios/heat3d_mpiio_con.c > CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.i

all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.s"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios/heat3d_mpiio_con.c -o CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.s

# Object files for target heat3d_mpiio_con
heat3d_mpiio_con_OBJECTS = \
"CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.o"

# External object files for target heat3d_mpiio_con
heat3d_mpiio_con_EXTERNAL_OBJECTS =

bin/heat3d_mpiio_con: all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/heat3d_mpiio_con.c.o
bin/heat3d_mpiio_con: all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/build.make
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/libpng-1.6.36/build/lib/libpng.a
bin/heat3d_mpiio_con: nmoments-anal/libnmoments_analysis.a
bin/heat3d_mpiio_con: all-transports/adios-all/libadaptor.a
bin/heat3d_mpiio_con: /usr/lib64/libpng.so
bin/heat3d_mpiio_con: /usr/lib64/libjpeg.so
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdecaf.a
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libmanala.a
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_datamodel.a
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbredala_transport_mpi.a
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libbca.a
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/code/decaf/install/lib/libdca.a
bin/heat3d_mpiio_con: /lustre/atlas/proj-shared/csc143/dhuang/sc15/libpng-1.6.36/build/lib/libpng.a
bin/heat3d_mpiio_con: all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable ../../../bin/heat3d_mpiio_con"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/heat3d_mpiio_con.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/build: bin/heat3d_mpiio_con

.PHONY : all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/build

all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/clean:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios && $(CMAKE_COMMAND) -P CMakeFiles/heat3d_mpiio_con.dir/cmake_clean.cmake
.PHONY : all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/clean

all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/depend:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : all-transports/adios-all/heat3d-adios/CMakeFiles/heat3d_mpiio_con.dir/depend

