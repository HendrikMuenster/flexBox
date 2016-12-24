set(INCLUDES_LOCAL)
list(APPEND INCLUDES_LOCAL
  flexDiagonalOperator.h
  flexGradientOperator.h
  flexIdentityOperator.h
  flexLinearOperator.h
  flexSuperpixelOperator.h
  flexZeroOperator.h
)

if(USE_CUDA)
  list(APPEND INCLUDES_LOCAL flexMatrixGPU.h)
else()
  list(APPEND INCLUDES_LOCAL flexMatrix.h)
endif(USE_CUDA)

#propagate to main cmake file
set(INCLUDES_OPERATOR ${INCLUDES_LOCAL} PARENT_SCOPE)
