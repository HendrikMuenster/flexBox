set(INCLUDES_LOCAL)
list(APPEND INCLUDES_LOCAL flexSolver.h)

if(USE_CUDA)
  list(APPEND INCLUDES_LOCAL flexSolverPrimalDualCuda.h)
else()
  list(APPEND INCLUDES_LOCAL flexSolverPrimalDual.h)
endif(USE_CUDA)

#propagate to main cmake file
set(INCLUDES_SOLVER ${INCLUDES_LOCAL} PARENT_SCOPE)
