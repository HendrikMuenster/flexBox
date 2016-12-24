#functions.{cpp | cu} is matlab specific (contains mex entry functions)
if(USE_CUDA)
  #set(SOURCES functions.cu)
else()
  matlab_add_mex (NAME flexBoxCPP
    SRC functions.cpp
    LINK_TO flexBoxLib ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${MAT_LIBRARY}
  )
endif(USE_CUDA)
