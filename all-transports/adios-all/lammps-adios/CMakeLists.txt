#if(LAMMPS_INCLUDE_DIR)
    #add_definitions("-DPRECISE")

    include_directories       (${LAMMPS_INCLUDE_DIR})

    add_executable            (lammps_adios_prod        lammps_adios_prod.cpp )
    add_executable            (lammps_adios_con     lammps_adios_con.cpp lammps_select.cpp)
    add_executable            (lammps_mpiio_con     lammps_mpiio_con.cpp lammps_select.cpp)

    add_executable            (lammps_native_con     lammps_native_con.cpp)

#    add_executable	      (lammps_sensei    lammpsBridge.cxx lammpsDataAdaptor.cxx 
#                               lammps_sensei_analysis.cpp lammps_sensei_test.cpp)
    add_executable            (lammps_sensei    lammps_sensei_analysis.cpp lammps_sensei_test.cpp)

    target_link_libraries     (lammps_adios_prod
                                ${LAMMPS_LIBRARY_STATIC}
				${ADIOS_LIB}
                                msd_analysis
                                adaptor m)

    target_link_libraries     (lammps_adios_con
				${ADIOS_LIB}
                                ${LAMMPS_LIBRARY_STATIC}
                                msd_analysis
                                adaptor m)
        # ${bredala_transports} bredala_datamodel manala decaf)

    target_link_libraries     (lammps_mpiio_con
				${ADIOS_LIB}
				${LAMMPS_LIBRARY_STATIC}
                                msd_analysis m
                                adaptor)

    target_link_libraries     (lammps_native_con
				${ADIOS_LIB}
				${LAMMPS_LIBRARY_STATIC}
                                msd_analysis m
                                adaptor)

    target_link_libraries     (lammps_sensei
                                ${LAMMPS_LIBRARY_STATIC}
				${ADIOS_LIB}
				sensei 
                                msd_analysis m) 


#endif()



