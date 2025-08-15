cmake_minimum_required(VERSION 3.22)
project(bipm_g2o LANGUAGES CXX)


# Ensure C++17 is used
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
 

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

#----------------------------------------------------------------------------------------
# Find required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(g2o 1.0.0 REQUIRED)




# Specify the path to g2o include directory and libraries
set(G2O_INCLUDE_DIRS "/usr/include/include" CACHE PATH "g2o include directory")
set(G2O_LIBRARY_DIRS "/usr/include/lib" CACHE PATH "g2o library directory")
set(G2O_LIBRARIES
    g2o::core
    g2o::stuff
    g2o::solver_cholmod
)

# OpenMP support
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DG2O_USE_OPENMP)
endif()


# Manually specify g2o include directories and libraries
include_directories(${G2O_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
link_directories(${G2O_LIBRARY_DIRS})

# Manually specify bipm_g2o include directories and libraries
include_directories(${G2O_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)  
link_directories(${G2O_LIBRARY_DIRS})


# add executable test
add_executable(example example/example.cpp)
target_link_libraries(example ${G2O_LIBRARIES} Eigen3::Eigen fmt::fmt)
#----------------------------------------------------------------------------------------


# Configuration option: Default to USE_BIPM, but allow switching to USE_AL
#option(USE_BIPM "Build with BIPM features (disable for AL )" ON)
option(USE_BIPM "Build with AL features (disable for BIPM )" ON)


if(USE_BIPM)
    message(STATUS "Building with USE_BIPM enabled")
    add_definitions(-DUSE_BIPM)
else()
    message(STATUS "Building with USE_AL enabled (USE_BIPM disabled)")
    add_definitions(-DUSE_AL)
endif()


