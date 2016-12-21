cmake_minimum_required (VERSION 3.1)
project (FlexBox)

set(CMAKE_CXX_STANDARD 11)


if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message (STATUS  "No build type specified. Default is 'Release'.")
        set (CMAKE_BUILD_TYPE Release CACHE STRING "The type of build" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif ()

#add_subdirectory(...)
#include_directories(...)

#TODO: build dirs
#TODO: debug/release
#TODO: python interface?
#TODO: status, errors, etc.
#TODO: CMAKE_BUILD_TYPE, CMAKE_CONFIGURATION_TYPES
#TODO: maybe add install target
#TODO: maybe add shared/static option
#TODO: MX, ENG, MAT - Libs?
#TODO: add support for older cmake versions?
#TODO: add option to build examples

option (USE_OPENMP "Use OpenMP" ON)
option (USE_CUDA "Use CUDA" OFF)
option (BUILD_MATLAB "build Matlab (mex) interface" ON)

#set (WRAPPERS_DIR "wrappers") #if f.ex python wrapper is added in the future


include_directories (.)

#includes and sources for all configs
#set(SOURCES main.cpp)
set(INCLUDES
  flexBox.h
  flexBoxData.h
  flexDiagonalOperator.h
  flexDualizedDataTerm.h
  flexDualizedOperator.h
  flexGradientOperator.h
  flexIdentityOperator.h
  flexLinearOperator.h
  flexProxList.h
  flexSolver.h
  flexSuperpixelOperator.h
  flexTermDual.h
  flexTermPrimal.h
  flexTildeMultiOperatorMultiDual.h
  flexZeroOperator.h
  tools.h
  )


#set additional build flags
set (DEBUG_FLAGS "-g -DVERBOSE")
set (RELEASE_FLAGS "-O3 -DNDEBUG")

set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")

#OpenMP
if (USE_OPENMP)
  find_package (OpenMP)

  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else ()
    message (STATUS "Couldn't find OpenMP. Building without...")
  endif (OPENMP_FOUND)

endif (USE_OPENMP)

#includes cmakes 3.7. version of FindMatlab which is not buggy
list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmakeModules)

#Matlab
if (BUILD_MATLAB)
  find_package(Matlab)

  if (Matlab_FOUND)

    #functions.{cpp | cu} is matlab specific (contains mex entry functions)
    if (USE_CUDA)
      #set(SOURCES functions.cu)
    else ()
      #set(SOURCES functions.cpp)
    endif (USE_CUDA)

  else ()
    message (STATUS "Couldn't find Matlab. Building without...")
  endif (Matlab_FOUND)

endif (BUILD_MATLAB)

if (USE_CUDA)
else ()
  #non cuda specific includes
  set(INCLUDES ${INCLUDES} flexBoxDataCPU.h flexMatrix.h flexSolverPrimalDual.h)
endif (USE_CUDA)

#MESSAGE ( STATUS "includes: " ${INCLUDES} )
#MESSAGE ( STATUS "SOURCES: " ${SOURCES} )

#build flexBox library
add_library(flexBoxLib SHARED ${INCLUDES} ${SOURCES})
target_include_directories(flexBoxLib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
set_target_properties(flexBoxLib PROPERTIES LINKER_LANGUAGE CXX)

#MESSAGE ( STATUS "Matlab_LIBRARIES: " ${Matlab_LIBRARIES} )
#MESSAGE ( STATUS "Matlab_MEX_LIBRARY: " ${Matlab_MEX_LIBRARY} )
#MESSAGE ( STATUS "Matlab_MX_LIBRARY: " ${Matlab_MX_LIBRARY} )
#MESSAGE ( STATUS "Matlab_ENG_LIBRARY: " ${Matlab_ENG_LIBRARY} )

if (Matlab_FOUND)
  matlab_add_mex (NAME flexBoxCPP2
    SRC functions.cpp
    LINK_TO ${Matlab_LIBRARIES} #flexBoxLib ${Matlab_LIBRARIES} why is flexBoxLib not needed?
  )
endif ()


#build
add_executable(flexBox main.cpp)
target_include_directories(flexBox PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(flexBox PUBLIC flexBoxLib)

#add_executable(FlexBox ${SOURCES} ${INCLUDES})
#target_link_libraries
