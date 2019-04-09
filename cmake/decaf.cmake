project                     (Decaf)
cmake_minimum_required      (VERSION 3.0)

option                      (transport_mpi      "Build Decaf with MPI transport layer"          ON)
option                      (transport_cci      "Build Decaf with CCI transport layer"          OFF)
option                      (transport_file     "Build Decaf with file transport layer"         OFF)
option                      (tess_dense         "Build tessellation density estimator example"  OFF)
option                      (build_bredala      "Build Bredala libraries and examples"          ON)
option                      (build_manala       "Build Manala libraries and examples"           ON)
option                      (build_decaf        "Build the Decaf workflow system"               ON)
option                      (build_tests        "Build the tests examples"                      ON)
set                         (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

if(build_manala)
    set                     (build_bredala  true)
endif(build_manala)
if(build_decaf)
    set                     (build_bredala  true)
    set                     (build_manala   true)
endif(build_decaf)

# OSX flags
if                          (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    add_definitions         (-DMAC_OSX)
    set                     (CMAKE_MACOSX_RPATH	  on)

    # --- following RPATH settings are for Sierra w/ Clang, hopefully they don't hurt other versions
    # ref: https://cmake.org/Wiki/CMake_RPATH_handling
    # use, i.e. don't skip, the full RPATH for the build tree
    set                     (CMAKE_SKIP_BUILD_RPATH            false)
    # when building, don't use the install RPATH already (but later on when installing)
    set                     (CMAKE_BUILD_WITH_INSTALL_RPATH    false)
    # set RPATH to install path
    set                     (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set                     (CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
    # the RPATH to be used when installing, but only if it's not a system directory
    list                    (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
                                                               "${CMAKE_INSTALL_PREFIX}/lib"
                                                               isSystemDir)
    if                      ("${isSystemDir}" STREQUAL         "-1")
      set                   (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    endif                   ()
endif                       (${CMAKE_SYSTEM_NAME} MATCHES      "Darwin")

#SET(CMAKE_C_COMPILER /opt/gcc/5.3.0/bin/gcc)
#SET(CMAKE_CXX_COMPILER /opt/gcc/5.3.0/bin/g++)
SET(CMAKE_C_COMPILER cc)
SET(CMAKE_CXX_COMPILER CC)

# C++11
#set                         (CMAKE_CXX_STANDARD        11)
set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11 -fPIC")
#set                         (CMAKE_CXX_STANDARD        11)
#set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# MPI
if                          (transport_mpi)
#find_package(MPI REQUIRED)
#if (MPI_FOUND)
	set (MPI_INCLUDE_PATH "$ENV{CRAY_MPICH2_DIR}/include")
	set (MPI_C_LIBRARYIES "$ENV{CRAY_MPICH2_DIR}/lib/libmpich.so")
	set (MPI_CXX_LIBRARYIES "$ENV{CRAY_MPICH2_DIR}/lib/libmpichcxx.so")
	message("MPI is found! 1 ${MPI_INCLUDE_PATH} 2 ${MPI_C_LIBRARYIES} 3 ${MPI_CXX_LIBRARYIES} 4 ${MPI_CXX_INCLUDE_PATH} 5 ${MPI_CXX_COMPILE_FLAGS} 6 ${MPI_CXX_LIBRARIES} 7 ${MPI_CXX_LINK_FLAGS} mpi path")
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
#else (MPI_FOUND)
#    message(SEND_ERROR "This application cannot compile without MPI")
#endif (MPI_FOUND)
  if                        (NOT bgq)
    set                     (transport_libraries    ${transport_libraries} ${MPI_C_LIBRARIES}  ${MPI_CXX_LIBRARIES})
  endif                     ()
  add_definitions           (-DTRANSPORT_MPI)
  set                       (TRANSPORT_MPI ON)
endif                       (transport_mpi)

#CCI
if                          (transport_cci)
  find_package              (MPI REQUIRED)
  find_package              (CCI REQUIRED)
  set                       (transport_libraries    ${transport_libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${CCI_LIBRARY})
  include_directories       (${CCI_INCLUDE_DIR})
  add_definitions           (-DTRANSPORT_CCI)
  set                       (TRANSPORT_CCI ON)
endif                       (transport_cci)

#FILE
#Should be used with the variable HDF5_PREFER_PARALLEL set to true
if                          (transport_file)
  find_package              (HDF5 REQUIRED)
  include_directories       (${HDF5_INCLUDE_DIR})
  set                       (transport_libraries    ${transport_libraries} ${HDF5_LIBRARIES})
  add_definitions           (-DTRANSPORT_FILE)
  set                       (TRANSPORT_FILE ON)
endif                       (transport_file)

# Boost
find_package                (Boost 1.59.0 COMPONENTS serialization python  REQUIRED)
message                     (STATUS "Boost libraries: " ${Boost_LIBRARIES})

# Set include directories
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${Boost_INCLUDE_DIRS}
                             ${CMAKE_CURRENT_BINARY_DIR}
                             ${CMAKE_CURRENT_SOURCE_DIR}
                             ${CMAKE_CURRENT_SOURCE_DIR}/include
                             ${MPI_INCLUDE_PATH})

# Set libraries
set                         (libraries
                             ${libraries}
                             ${transport_libraries}
                             ${CMAKE_DL_LIBS})

set (CMAKE_LINKER_FLAGS ${CMAKE_LINKER_FLAGS} "-fPIC -Wl,--export-dynamic -dynamic")

