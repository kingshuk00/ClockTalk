add_library(clock STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/clocks.c
)

target_include_directories(clock PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(clock PUBLIC paraver utils m)
