## setting  adios 
option(DEBUG_adaptor "debug option for dsadator and adios adaptor" OFF)
if(DEBUG_adaptor)
add_definitions(-DDEBUG)
endif(DEBUG_adaptor)
message("$ENV{DATASPACES_DIR}")
set(DS_DIR "$ENV{DATASPACES_DIR}")
set(ADIOS_DIR  "$ENV{ADIOS_DIR}")
set(FLEXPATH_DIR "$ENV{FLEXPATH_DIR}")
set(SENSEI_INSTALL "$ENV{SENSEI_INSTALL}")

#include_directories(${DS_DIR}/include ${ADIOS_DIR}/include ${FLEXPATH_DIR}/include /ccs/home/dhuang/sc15/enet-1.2.5/build/include /ccs/home/dhuang/sc15/bzip2-1.0.6-master/build/include $ENV{ADIOS_INC})
include_directories(${DS_DIR}/include ${ADIOS_DIR}/include ${FLEXPATH_DIR}/include $ENV{ADIOS_INC})
message("=====before find library ${DS_DIR}   ${DATASPACES_DIR}")
find_library(DSPACES  dspaces ${DS_DIR}/lib)
find_library(DSCOMMON dscommon ${DS_DIR}/lib) 
find_library(DART dart ${DS_DIR}/lib) 
find_library(ADIOS adios ${ADIOS_DIR}/lib)
#find_library(EVPATH evpath ${FLEXPATH_DIR}/lib)
#find_library(FFS ffs ${FLEXPATH_DIR}/lib)
#find_library(ATL atl ${FLEXPATH_DIR}/lib)
#find_library(DILL dill ${FLEXPATH_DIR}/lib)
#find_library(CERCS_ENV cercs_env ${FLEXPATH_DIR}/lib)
set(ENET "$ENV{MYPROJ}/sc15/enet-1.2.5/build/lib/libenet.a")
set(BZ2 "$ENV{MYPROJ}/sc15/bzip2-1.0.6-master/build/lib/libbz2.a")
set(PNG "$ENV{MYPROJ}/sc15/libpng-1.6.36/build/lib/libpng.a")
message("=== after find library PNG: ${PNG}")
#message("=== variables: ${ADIOS} ${DSPACES} ${DSCOMMON} ${DART} ${EVPATH} ${FFS} ${ATL} ${DILL} ${CERCS_ENV} ${ENET}")
#set(ADIOS_LIB ${ADIOS} $ENV{ADIOS_CLIB} $ENV{ADIOSREAD_CLIB} ${DSPACES} ${DSCOMMON} ${DART} ${EVPATH} ${FFS} ${ATL} ${DILL} ${CERCS_ENV} ${ENET} ${BZ2} $ENV{SZ_LIB})

#string(STRIP "$ENV{ADIOS_CLIB}" ADIOS_CLIB)
#string(STRIP "$ENV{ADIOSREAD_CLIB}" ADIOSREAD_CLIB)
set(ADIOS_CLIB "-L${ADIOS_DIR}/lib -ladios -ladiosread  -L${DS_DIR}/lib  -ldspaces -ldscommon -ldart -ldspaces -ldscommon -ldart -L${FLEXPATH_DIR}/lib -levpath -lffs -latl -ldill -lcercs_env -lenet   -lpmi -lugni -llustreapi -lz -lbz2  -lrca -lm -lrt $ENV{MXML_LIB} $ENV{SZIP_POST_LINK_OPTS} $ENV{SZ_POST_LINK_OPTS}")
set(ADIOS_LIB ${PNG}  ${JPEG_LIB_PATH} ${ADIOS_CLIB})
find_package(SENSEI REQUIRED)
set(SENSEI_LIB "${SENSEI_INSTALL}/lib/libsensei.a"  "${SENSEI_INSTALL}/lib/libsenseiCore.a" "${SENSEI_INSTALL}/lib/libpugixml.a" "${SENSEI_INSTALL}/lib/libtimer.a")
#set(ADIOS_LIB ${PNG} -L/sw/xk6/adios/1.13.0/cle5.2_pgi16.10.0/lib -ladios -L/opt/cray/pmi/default/lib64 -L/opt/cray/ugni/default/lib64 -L/sw/xk6/sz/1.4.11/cle5.2_pgi16.10.0/lib -L/sw/xk6/dataspaces/1.6.2/cle5.2_pgi16.10.0/lib -ldspaces -ldscommon -ldart -L/sw/xk6/dataspaces/1.6.2/cle5.2_pgi16.10.0/lib -ldspaces -ldscommon -ldart /sw/xk6/flexpath/1.13/cle5.2_pgi16.10.0/lib/libevpath.a /sw/xk6/flexpath/1.13/cle5.2_pgi16.10.0/lib/libffs.a /sw/xk6/flexpath/1.13/cle5.2_pgi16.10.0/lib/libatl.a /sw/xk6/flexpath/1.13/cle5.2_pgi16.10.0/lib/libdill.a /sw/xk6/flexpath/1.13/cle5.2_pgi16.10.0/lib/libenet.a -lpmi -lugni -llustreapi -lz -lbz2 -lSZ -lzlib)

#include_directories("${CMAKE_SOURCE_DIR}/src")    
 
string(STRIP "${ADIOS_LIB}" ADIOS_LIB)

message("=== before add adaptor: ${ADIOS_LIB}aaa")


add_library(adaptor adios_adaptor.c ds_adaptor.c)
target_link_libraries(adaptor pthread rt jpeg m ${ADIOS_LIB})

# add libm if using gcc
# Dan comment
#target_link_libraries(adaptor  pthread rt rdmacm ibverbs)
#target_link_libraries(adaptor pthread rt jpeg m)

message("=== after add adaptor")
include_directories(${PROJECT_SOURCE_DIR}/all-transports/adios-all)
message("=== after include directory")
#add_subdirectory(lbm-adios)

#if(LAMMPS_INCLUDE_DIR)
add_subdirectory(lammps-adios)
#endif()

add_subdirectory(laplace-adios)
