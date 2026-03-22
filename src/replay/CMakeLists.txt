add_library(replay STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/replay.c
)

target_include_directories(replay PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(replay PUBLIC monitoring clock paraver trace utils m)
