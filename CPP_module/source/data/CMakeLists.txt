set(INCLUDES_LOCAL flexBoxData.h)

if(USE_CUDA)
  list(APPEND INCLUDES_LOCAL flexBoxDataGPU.h)
else()
  list(APPEND INCLUDES_LOCAL flexBoxDataCPU.h)
endif(USE_CUDA)

#propagate to main cmake file
set(INCLUDES_DATA ${INCLUDES_LOCAL} PARENT_SCOPE)
