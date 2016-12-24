set(INCLUDES_LOCAL)
list(APPEND INCLUDES_LOCAL
  flexProx.h
  flexProxDualBoxConstraint.h
  flexProxDualDataKL.h
  flexProxDualDataL1.h
  flexProxDualDataL2.h
  flexProxDualFrobenius.h
  flexProxDualHuber.h
  flexProxDualL1Ansio.h
  flexProxDualL1Iso.h
  flexProxDualL2.h
)

#propagate to main cmake file
set(INCLUDES_PROX ${INCLUDES_LOCAL} PARENT_SCOPE)
