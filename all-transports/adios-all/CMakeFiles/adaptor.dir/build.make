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
include all-transports/adios-all/CMakeFiles/adaptor.dir/depend.make

# Include the progress variables for this target.
include all-transports/adios-all/CMakeFiles/adaptor.dir/progress.make

# Include the compile flags for this target's objects.
include all-transports/adios-all/CMakeFiles/adaptor.dir/flags.make

all-transports/adios-all/CMakeFiles/adaptor.dir/adios_adaptor.c.o: all-transports/adios-all/CMakeFiles/adaptor.dir/flags.make
all-transports/adios-all/CMakeFiles/adaptor.dir/adios_adaptor.c.o: all-transports/adios-all/adios_adaptor.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object all-transports/adios-all/CMakeFiles/adaptor.dir/adios_adaptor.c.o"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/adaptor.dir/adios_adaptor.c.o   -c /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/adios_adaptor.c

all-transports/adios-all/CMakeFiles/adaptor.dir/adios_adaptor.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/adaptor.dir/adios_adaptor.c.i"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/adios_adaptor.c > CMakeFiles/adaptor.dir/adios_adaptor.c.i

all-transports/adios-all/CMakeFiles/adaptor.dir/adios_adaptor.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/adaptor.dir/adios_adaptor.c.s"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/adios_adaptor.c -o CMakeFiles/adaptor.dir/adios_adaptor.c.s

all-transports/adios-all/CMakeFiles/adaptor.dir/ds_adaptor.c.o: all-transports/adios-all/CMakeFiles/adaptor.dir/flags.make
all-transports/adios-all/CMakeFiles/adaptor.dir/ds_adaptor.c.o: all-transports/adios-all/ds_adaptor.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object all-transports/adios-all/CMakeFiles/adaptor.dir/ds_adaptor.c.o"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/adaptor.dir/ds_adaptor.c.o   -c /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/ds_adaptor.c

all-transports/adios-all/CMakeFiles/adaptor.dir/ds_adaptor.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/adaptor.dir/ds_adaptor.c.i"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/ds_adaptor.c > CMakeFiles/adaptor.dir/ds_adaptor.c.i

all-transports/adios-all/CMakeFiles/adaptor.dir/ds_adaptor.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/adaptor.dir/ds_adaptor.c.s"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && /opt/cray/craype/2.5.13/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/ds_adaptor.c -o CMakeFiles/adaptor.dir/ds_adaptor.c.s

# Object files for target adaptor
adaptor_OBJECTS = \
"CMakeFiles/adaptor.dir/adios_adaptor.c.o" \
"CMakeFiles/adaptor.dir/ds_adaptor.c.o"

# External object files for target adaptor
adaptor_EXTERNAL_OBJECTS =

all-transports/adios-all/libadaptor.a: all-transports/adios-all/CMakeFiles/adaptor.dir/adios_adaptor.c.o
all-transports/adios-all/libadaptor.a: all-transports/adios-all/CMakeFiles/adaptor.dir/ds_adaptor.c.o
all-transports/adios-all/libadaptor.a: all-transports/adios-all/CMakeFiles/adaptor.dir/build.make
all-transports/adios-all/libadaptor.a: all-transports/adios-all/CMakeFiles/adaptor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library libadaptor.a"
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && $(CMAKE_COMMAND) -P CMakeFiles/adaptor.dir/cmake_clean_target.cmake
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/adaptor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
all-transports/adios-all/CMakeFiles/adaptor.dir/build: all-transports/adios-all/libadaptor.a

.PHONY : all-transports/adios-all/CMakeFiles/adaptor.dir/build

all-transports/adios-all/CMakeFiles/adaptor.dir/clean:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all && $(CMAKE_COMMAND) -P CMakeFiles/adaptor.dir/cmake_clean.cmake
.PHONY : all-transports/adios-all/CMakeFiles/adaptor.dir/clean

all-transports/adios-all/CMakeFiles/adaptor.dir/depend:
	cd /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/all-transports/adios-all/CMakeFiles/adaptor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : all-transports/adios-all/CMakeFiles/adaptor.dir/depend

