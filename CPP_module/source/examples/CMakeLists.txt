#TODO: install target vs configurable files for data files

configure_file(
  sampleFiles.h.in
  sampleFiles.h
)

if(CUDA_FOUND)

  set(SAMPLE_SRCS_CUDA
    exampleROF.cu
  )

  foreach(SAMPLE_FILE ${SAMPLE_SRCS_CUDA})

    get_filename_component(FILE_NAME ${SAMPLE_FILE} NAME_WE)
    include_directories(${INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    cuda_add_executable(${FILE_NAME} ${SAMPLE_FILE})

    install(TARGETS ${FILE_NAME} DESTINATION ${CMAKE_BINARY_DIR}/bin)

  endforeach()

else()

  set(SAMPLE_SRCS_CPP
    exampleROF.cpp
  )
  foreach(SAMPLE_FILE ${SAMPLE_SRCS_CPP})

    get_filename_component(FILE_NAME ${SAMPLE_FILE} NAME_WE)
    add_executable(${FILE_NAME} ${SAMPLE_FILE})
    target_link_libraries(${FILE_NAME} PUBLIC flexBoxLib)
    target_include_directories(${FILE_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
    target_include_directories(${FILE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    install(TARGETS ${FILE_NAME} DESTINATION ${CMAKE_BINARY_DIR}/bin)

  endforeach()

endif(CUDA_FOUND)
