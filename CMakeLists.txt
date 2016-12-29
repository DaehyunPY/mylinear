# CMake project file for mylib
cmake_minimum_required (VERSION 2.6.4)
project (mylib)
enable_language (Fortran)

set (LIB_NAME mylinear)
set (VERSION 150326)





include_directories (/usr/include)
include_directories (/usr/local/include)
include_directories ($ENV{HOME}/include)
link_directories (/usr/lib)
link_directories (/usr/local/lib)
link_directories ($ENV{HOME}/lib)

ADD_LIBRARY (${LIB_NAME}        SHARED ${LIB_NAME}.f90)
ADD_LIBRARY (${LIB_NAME}-static STATIC ${LIB_NAME}.f90)
target_link_libraries (${LIB_NAME}        -llapack)
target_link_libraries (${LIB_NAME}-static -llapack)
set_target_properties ( ${LIB_NAME}
    PROPERTIES
    OUTPUT_NAME ${LIB_NAME}
    PREFIX lib
    VERSION ${VERSION}
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/../include )
set_target_properties ( ${LIB_NAME}-static
    PROPERTIES
    OUTPUT_NAME ${LIB_NAME}
    PREFIX lib
    VERSION ${VERSION}
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/../include )





### make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
            "Choose the type of build, options are: None Debug Release."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

### FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    ### gfortran
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    ### ifort (untested)
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
    ### g77
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
