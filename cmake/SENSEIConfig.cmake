set(SENSEI_VERSION "v0.0.0")
set(SENSEI_VERSION_MAJOR "0")
set(SENSEI_VERSION_MINOR "0")
set(SENSEI_VERSION_PATCH "0")
set(SENSEI_VERSION_DEVEL "")

if (NOT SENSEI_DIR)
  if (OFF)
    set(SENSEI_DIR "/lustre/atlas/proj-shared/csc143/zq53/sensei/sensei/build")
  else()
    set(SENSEI_DIR "/lustre/atlas/proj-shared/csc143/zq53/sensei/sensei/sensei-install")
  endif()
endif()
list(APPEND CMAKE_MODULE_PATH "${SENSEI_DIR}")

set(SENSEI_LIB_TYPE STATIC)
if (OFF)
  set(SENSEI_LIB_TYPE SHARED)
endif()

if (OFF)
  set(BUILD_SHARED_LIBS OFF FORCE)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(LINK_SEARCH_START_STATIC TRUE)
  set(LINK_SEARCH_END_STATIC TRUE)
endif()

set(ENABLE_SENSEI ON)
set(ENABLE_PYTHON OFF)
set(ENABLE_CATALYST OFF)
set(ENABLE_CATALYST_PYTHON OFF)
set(ENABLE_LIBSIM OFF)
set(ENABLE_ADIOS OFF)
set(ENABLE_CONDUIT OFF)
set(ENABLE_VTK_GENERIC_ARRAYS OFF)
set(ENABLE_VTK_ACCELERATORS OFF)
set(ENABLE_VTK_MPI OFF)
set(ENABLE_VTK_IO OFF)
set(ENABLE_VTKM OFF)
set(ENABLE_VTKM_RENDERING OFF)

if (NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "-fpermissive -fPIC -std=c++11 -Wall -Wextra -O3 -march=native -mtune=native"
    CACHE STRING "sensei build defaults"
  FORCE)
endif()

if (ENABLE_CATALYST)
  if (NOT ParaView_DIR)
    set(ParaView_DIR "")
    find_package(ParaView REQUIRED)
  endif()
else()
  if (NOT VTK_DIR)
    set(VTK_DIR "/sw/xk6/vtk/7.0.0/sles11.3_gnu4.9.0/lib/cmake/vtk-7.0")
  endif()
  find_package(VTK REQUIRED)
endif()

include(thread)
include(sMPI)
include(sVTK)
include(timer)
include(pugixml)
include(timer)
include(sDIY)
if (ENABLE_VTKM)
  include(sVTKm)
endif()
if (ENABLE_LIBSIM)
  include(sLibsim)
endif()
if (ENABLE_ADIOS)
  include(sADIOS)
endif()
if (ENABLE_CONDUIT)
  include(sConduit)
endif()
include(senseiCore)
if (ENABLE_PYTHON)
  include(sPython)
  include(_PythonAnalysis)
endif()
include(sensei)