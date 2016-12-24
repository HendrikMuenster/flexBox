cmake_minimum_required(VERSION 3.1)
project(FlexBox)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message (STATUS  "No build type specified. Default is 'Release'.")
        set (CMAKE_BUILD_TYPE Release CACHE STRING "The type of build" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

#TODO: build dirs
#TODO: python interface?
#TODO: maybe add shared lib later
#TODO: add support for older cmake versions?
#TODO: add option to build examples

option(USE_OPENMP "Use OpenMP" OFF)
option(USE_CUDA "Use CUDA" OFF)
option(BUILD_MATLAB "Build Matlab (mex) interface" ON)
option(BUILD_EXAMPLES "Build C++ examples" OFF)

set (WRAPPERS_DIR "wrappers")

#set additional build flags
set(DEBUG_FLAGS "-g -DVERBOSE")
set(RELEASE_FLAGS "-O3 -DNDEBUG")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")

macro( set_subdir_sources )
    get_directory_property( _PARENTDIR PARENT_DIRECTORY )
    file( RELATIVE_PATH _SUBDIR ${_PARENTDIR} ${CMAKE_CURRENT_SOURCE_DIR} )
    foreach( SOURCE_FILE ${ARGN} )
        set( SUBDIR_SOURCES ${SUBDIR_SOURCES} ${_SUBDIR}/${SOURCE_FILE} )
    endforeach()
    set( SUBDIR_SOURCES ${SUBDIR_SOURCES} PARENT_SCOPE )
endmacro()

add_subdirectory(data)
add_subdirectory(operator)
add_subdirectory(prox)
add_subdirectory(solver)
add_subdirectory(term)
set(INCLUDES flexBox.h tools.h)

#OpenMP
if(USE_OPENMP)
  find_package(OpenMP)

  if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    message (STATUS "Couldn't find OpenMP. Building without...")
  endif(OPENMP_FOUND)

endif(USE_OPENMP)

#includes FindMatlab (from cmake v3.7)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

#Matlab
if(BUILD_MATLAB)
  find_package(Matlab REQUIRED MX_LIBRARY MAT_LIBRARY)

  if(NOT Matlab_FOUND)
    message (STATUS "Couldn't find Matlab or one of its required components (MX and/or MAT). Building without...")
  endif()

endif(BUILD_MATLAB)

#CUDA
if(USE_CUDA)
  find_package(CUDA)

  if(NOT CUDA_FOUND)
    message (STATUS "Couldn't find CUDA. Building without...")
  endif()

endif(USE_CUDA)

set(INCLUDE_DIRS)
list(APPEND INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/data
  ${CMAKE_CURRENT_SOURCE_DIR}/operator
  ${CMAKE_CURRENT_SOURCE_DIR}/prox
  ${CMAKE_CURRENT_SOURCE_DIR}/solver
  ${CMAKE_CURRENT_SOURCE_DIR}/term
)
#flexBox template library
add_library(flexBoxLib INTERFACE)
target_include_directories(flexBoxLib INTERFACE ${INCLUDE_DIRS})
#set_target_properties(flexBoxLib PROPERTIES LINKER_LANGUAGE CXX)
#set_target_properties(flexBoxLib PROPERTIES OUTPUT_NAME flexBox)

#build main
add_executable(flexBoxMain main.cpp)
#target_include_directories(flexBox PUBLIC INCLUDE_DIRS)
target_link_libraries(flexBoxMain PUBLIC flexBoxLib)


add_custom_target(flexBoxIDE)
target_sources()


#build Matlab
if(Matlab_FOUND)
  add_subdirectory(${WRAPPERS_DIR}/Matlab)
endif(Matlab_FOUND)
