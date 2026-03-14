add_library(monitoring STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/monitoring.c
)

target_include_directories(monitoring PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(monitoring PUBLIC trace m)
