#if(LAMMPS_INCLUDE_DIR)
    #add_definitions("-DPRECISE")

    include_directories       (.)

    add_executable            (heat3d_adios_prod        heat3d_adios_prod.c  explUtilPar.c updateBound.c readParam.c)
    add_executable            (heat3d_adios_con     heat3d_adios_con.c)
    add_executable            (heat3d_mpiio_con     heat3d_mpiio_con.c)
    add_executable            (heat3d_native_con     heat3d_native_con.c)
    add_library				(heat3d STATIC explicitPar.c explUtilPar.c updateBound.c readParam.c)
#    add_executable            (lammps_mpiio_con     lammps_mpiio_con.cpp lammps_select.cpp)
    target_link_libraries     (heat3d_adios_prod
				${ADIOS_LIB}
								nmoments_analysis
                                adaptor m)

    target_link_libraries     (heat3d_adios_con
				${ADIOS_LIB}
								nmoments_analysis
                                adaptor m)

    target_link_libraries     (heat3d
                ${ADIOS_LIB})

    target_link_libraries     (heat3d_mpiio_con
				${ADIOS_LIB}
								nmoments_analysis
                                adaptor m)
    target_link_libraries     (heat3d_native_con
				${ADIOS_LIB}
								nmoments_analysis
                                adaptor m)
foreach(item heat3d)

    #target_link_libraries(${item} ${MPI_C_LIBRARIES})

    if(MPI_COMPILE_FLAGS)
        set_target_properties(${item} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
    endif()

    if(MPI_LINK_FLAGS)
        set_target_properties(${item} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
    endif()
endforeach(item)

#	find_library(HEAT3D_LIB_PATH NAMES heat3d PATHS .) 
#set(common_library ${common_library} ${HEAT3D_LIB_PATH})
        # ${bredala_transports} bredala_datamodel manala decaf)

   # target_link_libraries     (lammps_mpiio_con
	#			${ADIOS_LIB}
	#			${LAMMPS_LIBRARY_STATIC}
     #                           msd_analysis m
      #                          adaptor
       #                         ${common_library})




#endif()



