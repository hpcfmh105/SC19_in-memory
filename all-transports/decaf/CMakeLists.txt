
message("Building vector example")
set(transport_mpi on CACHE BOOL "enable mpi transport for decaf")
#set(transport_file off CACHE BOOL "enable hdf5 file transport for decaf")
#set(CCI_PREFIX "${FWORK}/software/install" CACHE PATH "cci path")
include(${CMAKE_SOURCE_DIR}/cmake/decaf.cmake)
#include_directories(${CMAKE_SOURCE_DIR}/heat3d)
include_directories(${DECAF_INCLUDE_DIR})
include_directories(${HEAT3D_INCLUDE_DIR})
# heat3d

add_executable            (laplace_decaf     laplace_decaf.cpp  readParam.hpp)

add_executable            (laplace_in     laplace_in.cpp  readParam.hpp)
add_executable            (laplace_out     laplace_out.cpp  readParam.hpp)
add_library               (mod_laplace_decaf_shared MODULE laplace_decaf.cpp  readParam.hpp)
target_link_libraries     (laplace_decaf
                            nmoments_analysis
                            ${libraries}
                            ${common_library})

target_link_libraries     (laplace_in
                            nmoments_analysis
                            ${libraries}
                            ${common_library})
target_link_libraries     (laplace_out
                            nmoments_analysis
                            ${libraries}
                            ${common_library})
target_link_libraries   (mod_laplace_decaf_shared
						nmoments_analysis
                        ${libraries}
						${common_library})
 

   #${bredala_transports} bredala_datamodel manala decaf)

# .SO file extension on Linux/Mac OS
#set_target_properties(mod_heat3d_decaf PROPERTIES SUFFIX ".so")

# Don't add a 'lib' prefix to the shared library
#set_target_properties(mod_heat3d_decaf PROPERTIES PREFIX "")


#lammps
#find_path                   (LAMMPS_INCLUDE_DIR		lammps.h)
#find_library                (LAMMPS_LIBRARY NAMES	lammps)


if(LAMMPS_INCLUDE_DIR)
    include_directories       (${LAMMPS_INCLUDE_DIR})
    add_executable            (lammps_decaf            lammps_decaf.cpp)
    add_executable            (lammps_in            lammps_in.cpp)
    add_executable            (lammps_out            lammps_out.cpp)
    add_library               (mod_lammps_decaf_static STATIC lammps_decaf.cpp)
    add_library               (mod_lammps_decaf_shared SHARED lammps_decaf.cpp)
    target_link_libraries     (lammps_decaf
        ${LAMMPS_LIBRARY_STATIC}
                                msd_analysis
                                ${libraries}
                                ${common_library})
        # ${bredala_transports} bredala_datamodel manala decaf)
    target_link_libraries     (lammps_in
        ${LAMMPS_LIBRARY_STATIC}
                                msd_analysis
                                ${libraries}
                                ${common_library})
    target_link_libraries     (lammps_out
        ${LAMMPS_LIBRARY_STATIC}
                                msd_analysis
                                ${libraries}
                                ${common_library})

    target_link_libraries   (mod_lammps_decaf_static
                            ${libraries})

    target_link_libraries   (mod_lammps_decaf_shared
							msd_analysis_shared
							${LAMMPS_LIBRARY_SHARED}
							${common_library}
                            ${libraries})

  #set_target_properties(mod_lammps_decaf_shared PROPERTIES
   #     LINK_FLAGS "-shared -Wl,-soname,libmod_lammps_decaf_shared.so"
    #)
  # target_link_libraries   (mod_lammps_decaf_shared 
#							mod_lammps_decaf_static
 #                           ${libraries})

# Dan Modify,after run cmake ../CMakeList.txt, and make. Have to run the following code to generate loadable shared lib.
# /opt/cray/craype/2.5.13/bin/CC -fPIC -g -fopenmp -O3 -DDEBUG -DUSE_SAME_LOCK  -shared -Wl,-soname,libmod_lammps_decaf_shared.so  -o ../../lib/libmod_lammps_decaf_shared.so  libmod_lammps_decaf_static.a

        #${bredala_transports} bredala_datamodel manala decaf)

    # .SO file extension on Linux/Mac OS
    #set_target_properties(mod_lammps_decaf_shared PROPERTIES SUFFIX ".so")

    # Don't add a 'lib' prefix to the shared library
 #   set_target_properties(mod_lammps_decaf PROPERTIES PREFIX "")

endif()





