add_library(trace STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/trace.c
)

target_include_directories(trace PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(trace PUBLIC utils m)
