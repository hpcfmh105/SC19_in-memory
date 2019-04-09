# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/autofs/nccs-svm1_sw/titan/.swci/0-login/opt/spack/20180315/linux-suse_linux11-x86_64/gcc-5.3.0/cmake-3.11.3-f2xxt5vasvlzrefneolaxu4rgko7cjnb/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/autofs/nccs-svm1_sw/titan/.swci/0-login/opt/spack/20180315/linux-suse_linux11-x86_64/gcc-5.3.0/cmake-3.11.3-f2xxt5vasvlzrefneolaxu4rgko7cjnb/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named nmoments_analysis

# Build rule for target.
nmoments_analysis: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 nmoments_analysis
.PHONY : nmoments_analysis

# fast build rule for target.
nmoments_analysis/fast:
	$(MAKE) -f nmoments-anal/CMakeFiles/nmoments_analysis.dir/build.make nmoments-anal/CMakeFiles/nmoments_analysis.dir/build
.PHONY : nmoments_analysis/fast

#=============================================================================
# Target rules for targets named msd_analysis

# Build rule for target.
msd_analysis: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 msd_analysis
.PHONY : msd_analysis

# fast build rule for target.
msd_analysis/fast:
	$(MAKE) -f msd-anal/CMakeFiles/msd_analysis.dir/build.make msd-anal/CMakeFiles/msd_analysis.dir/build
.PHONY : msd_analysis/fast

#=============================================================================
# Target rules for targets named heat3d_decaf

# Build rule for target.
heat3d_decaf: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 heat3d_decaf
.PHONY : heat3d_decaf

# fast build rule for target.
heat3d_decaf/fast:
	$(MAKE) -f all-transports/decaf/CMakeFiles/heat3d_decaf.dir/build.make all-transports/decaf/CMakeFiles/heat3d_decaf.dir/build
.PHONY : heat3d_decaf/fast

#=============================================================================
# Target rules for targets named lammps_decaf

# Build rule for target.
lammps_decaf: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 lammps_decaf
.PHONY : lammps_decaf

# fast build rule for target.
lammps_decaf/fast:
	$(MAKE) -f all-transports/decaf/CMakeFiles/lammps_decaf.dir/build.make all-transports/decaf/CMakeFiles/lammps_decaf.dir/build
.PHONY : lammps_decaf/fast

#=============================================================================
# Target rules for targets named mod_heat3d_decaf

# Build rule for target.
mod_heat3d_decaf: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mod_heat3d_decaf
.PHONY : mod_heat3d_decaf

# fast build rule for target.
mod_heat3d_decaf/fast:
	$(MAKE) -f all-transports/decaf/CMakeFiles/mod_heat3d_decaf.dir/build.make all-transports/decaf/CMakeFiles/mod_heat3d_decaf.dir/build
.PHONY : mod_heat3d_decaf/fast

#=============================================================================
# Target rules for targets named mod_lammps_decaf_static

# Build rule for target.
mod_lammps_decaf_static: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mod_lammps_decaf_static
.PHONY : mod_lammps_decaf_static

# fast build rule for target.
mod_lammps_decaf_static/fast:
	$(MAKE) -f all-transports/decaf/CMakeFiles/mod_lammps_decaf_static.dir/build.make all-transports/decaf/CMakeFiles/mod_lammps_decaf_static.dir/build
.PHONY : mod_lammps_decaf_static/fast

#=============================================================================
# Target rules for targets named test_lammps_mpi

# Build rule for target.
test_lammps_mpi: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_lammps_mpi
.PHONY : test_lammps_mpi

# fast build rule for target.
test_lammps_mpi/fast:
	$(MAKE) -f tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/build.make tests/test-lammps/CMakeFiles/test_lammps_mpi.dir/build
.PHONY : test_lammps_mpi/fast

#=============================================================================
# Target rules for targets named test_lammps

# Build rule for target.
test_lammps: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_lammps
.PHONY : test_lammps

# fast build rule for target.
test_lammps/fast:
	$(MAKE) -f tests/test-lammps/CMakeFiles/test_lammps.dir/build.make tests/test-lammps/CMakeFiles/test_lammps.dir/build
.PHONY : test_lammps/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... nmoments_analysis"
	@echo "... msd_analysis"
	@echo "... heat3d_decaf"
	@echo "... lammps_decaf"
	@echo "... mod_heat3d_decaf"
	@echo "... mod_lammps_decaf_static"
	@echo "... test_lammps_mpi"
	@echo "... test_lammps"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
