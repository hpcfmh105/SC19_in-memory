
IF($ENV{CRAYPE_LINK_TYPE} MATCHES "dynamic")
add_subdirectory(decaf)
ELSE()
add_subdirectory(adios-all)
ENDIF($ENV{CRAYPE_LINK_TYPE} MATCHES "dynamic")
