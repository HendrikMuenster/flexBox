cmake_minimum_required(VERSION 3.1)
project(FlexBox)

#TODO: add support for older cmake versions?
#TODO: add option to build examples
#TODO: largeArrayDims??

##########################
#     FILES/INCUDES      #
##########################

#includes
set(INCLUDE_DIRS
  .
  data
  operator
  prox
  solver
  term
)

#data
set(SOURCES_IDE_DATA
  data/flexBoxData.h
  data/flexBoxDataCPU.h
  data/flexBoxDataGPU.h
)

#operator
set(SOURCES_IDE_OP
  operator/flexDiagonalOperator.h
  operator/flexGradientOperator.h
  operator/flexIdentityOperator.h
  operator/flexLinearOperator.h
  operator/flexMatrix.h
  operator/flexMatrixGPU.h
  operator/flexSuperPixelOperator.h
  operator/flexZeroOperator.h
)

#prox
set(SOURCES_IDE_PROX
  prox/flexProx.h
  prox/flexProxDualBoxConstraint.h
  prox/flexProxDualDataKL.h
  prox/flexProxDualDataL1.h
  prox/flexProxDualDataL2.h
  prox/flexProxDualFrobenius.h
  prox/flexProxDualHuber.h
  prox/flexProxDualL1Aniso.h
  prox/flexProxDualL1Iso.h
  prox/flexProxDualL2.h
)

#solver
set(SOURCES_IDE_SOLV
  solver/flexSolver.h
  solver/flexSolverPrimalDual.h
  solver/flexSolverPrimalDualCuda.h
)

#term
set(SOURCES_IDE_TERM
  term/flexTermDual.h
  term/flexTermPrimal.h
)

#matlab
set(SOURCES_IDE_MATLAB
  wrappers/Matlab/functions.cpp
  wrappers/Matlab/functions.cu
)

#group for visual studio
source_group("data" FILES ${SOURCES_IDE_DATA})
source_group("operator" FILES ${SOURCES_IDE_OP})
source_group("prox" FILES ${SOURCES_IDE_PROX})
source_group("solver" FILES ${SOURCES_IDE_SOLV})
source_group("term" FILES ${SOURCES_IDE_TERM})
source_group("matlab" FILES ${SOURCES_IDE_MATLAB})

#combine
set(SOURCES_IDE ${SOURCES_IDE_DATA} ${SOURCES_IDE_OP} ${SOURCES_IDE_PROX} ${SOURCES_IDE_SOLV} ${SOURCES_IDE_TERM} ${SOURCES_IDE_MATLAB})

##########################
#   END FILES/INCUDES    #
##########################

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message (STATUS  "No build type specified. Default is 'Release'.")
        set (CMAKE_BUILD_TYPE Release CACHE STRING "The type of build" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()


option(USE_OPENMP "Use OpenMP" OFF)
option(USE_CUDA "Use CUDA" OFF)
option(BUILD_MATLAB "Build Matlab (mex) interface" ON)
option(BUILD_EXAMPLES "Build C++/CUDA examples" OFF)

set (WRAPPERS_DIR "wrappers")

#set additional build flags
set(DEBUG_FLAGS "-g -DVERBOSE")
set(RELEASE_FLAGS "-O3 -DNDEBUG")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")

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
  find_package(Matlab COMPONENTS MX_LIBRARY MAT_LIBRARY)

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


#flexBox template library
add_library(flexBoxLib INTERFACE)
target_include_directories(flexBoxLib INTERFACE ${INCLUDE_DIRS})

#build main
add_executable(Example main.cpp)
target_link_libraries(Example PUBLIC flexBoxLib)


add_custom_target(IDE SOURCES ${SOURCES_IDE})


#add Matlab target
if(Matlab_FOUND)
  add_subdirectory(${WRAPPERS_DIR}/matlab)
endif(Matlab_FOUND)
