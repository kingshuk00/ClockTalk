#
# Copyright (c) 2026      Kingshuk Haldar.
#                         All rights reserved.
#
# Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
#

add_subdirectory(trace)
add_subdirectory(paraver)
add_subdirectory(clocks)
add_subdirectory(monitoring)

file(GLOB REPLAY_ROOT_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
)
list(FILTER REPLAY_ROOT_SOURCES EXCLUDE REGEX ".*/replay\\.c$")

add_library(replay STATIC ${REPLAY_ROOT_SOURCES})

target_include_directories(replay PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(replay PUBLIC
  paraver
  clock
  monitoring
  trace
  m
)


# main function
add_executable(clocktalk
  ${CMAKE_CURRENT_SOURCE_DIR}/replay.c
)

target_include_directories(clocktalk PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Static link order matters; dependencies should appear to the right.
target_link_libraries(clocktalk PRIVATE
  replay
  m
)
