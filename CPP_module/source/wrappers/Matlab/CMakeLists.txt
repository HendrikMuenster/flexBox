#functions.{cpp | cu} is matlab specific (contains mex entry functions)
if(USE_CUDA)
  #set(SOURCES functions.cu)
  cuda_add_executable(Matlab_CUDA functions.cu)
  target_link_libraries(Matlab_CUDA flexBoxLib ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${MAT_LIBRARY})

  matlab_add_mex(NAME)
else()
  matlab_add_mex (NAME Matlab
    SRC functions.cpp
    LINK_TO flexBoxLib ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${MAT_LIBRARY}
    OUTPUT_NAME flexBoxCPP
  )
endif(USE_CUDA)
