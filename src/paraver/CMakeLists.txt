add_library(paraver STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/paraver.c
)

target_include_directories(paraver PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(paraver PUBLIC trace m)
