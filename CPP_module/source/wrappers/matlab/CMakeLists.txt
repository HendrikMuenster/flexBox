#functions.{cpp | cu} is matlab specific (contains mex entry functions)
if(USE_CUDA)
  #must use include_directories instead of target_include_directories
  #known issue in cmakes CUDA-Module
  include_directories(${Matlab_INCLUDE_DIRS})
  include_directories(${INCLUDE_DIRS})
  set(CUDA_NVCC_FLAGS --std=c++11 -Xcompiler -fPIC -arch=compute_50 -code=sm_50)
  CUDA_COMPILE(cuda_compiled functions.cu SHARED)
  matlab_add_mex(NAME Matlab_CUDA
        SRC ${cuda_compiled}
        LINK_TO flexBoxLib ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${Matlab_MAT_LIBRARY} ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_CUFFT_LIBRARIES}
        OUTPUT_NAME flexBoxCPPCUDA
  )
  set_target_properties(Matlab_CUDA PROPERTIES LINKER_LANGUAGE CXX)

else()
  matlab_add_mex (NAME Matlab
    SRC functions.cpp
    LINK_TO flexBoxLib ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${MAT_LIBRARY}
    OUTPUT_NAME flexBoxCPP
  )
endif(USE_CUDA)
