set(INCLUDES_LOCAL)
list(APPEND INCLUDES_LOCAL flexTermDual.h flexTermPrimal.h)

#propagate to main cmake file
set(INCLUDES_TERM ${INCLUDES_LOCAL} PARENT_SCOPE)
