#
# Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
#
# Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
#

add_subdirectory(utils)
add_subdirectory(trace)
add_subdirectory(paraver)
add_subdirectory(clocks)
add_subdirectory(monitoring)
add_subdirectory(replay)

# link-time optimisation
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# main function
add_executable(clocktalk
  ${CMAKE_CURRENT_SOURCE_DIR}/opt_parser.c
  ${CMAKE_CURRENT_SOURCE_DIR}/clocktalk.c
)

target_include_directories(clocktalk PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Static link order matters; dependencies should appear to the right.
target_link_libraries(clocktalk PRIVATE
  replay
  paraver
  clock
  monitoring
  trace
  utils
  m
)
