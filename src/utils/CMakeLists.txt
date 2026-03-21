add_library(utils STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/utils.c
)

target_include_directories(utils PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(utils PUBLIC m)
