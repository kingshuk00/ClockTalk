#
# Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
#
# Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
#

cmake_minimum_required(VERSION 3.15)
project(ClockTalk VERSION 1.0.0 LANGUAGES C)

# Language options
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release/Debug)" FORCE)
endif()

# Compiler flags based on Release vs Debug
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG")   # optimize O3, define NDEBUG
set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG} -g3 -O0")          # debug info, no optimization

# Global warnings
add_compile_options(-Wall -Werror)   # all warnings are errors
# add_compile_options(-Wextra)         # extra warnings

option(CLOCKTALK_ENABLE_ASAN  "Enable AddressSanitizer" OFF)
option(CLOCKTALK_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
if(CLOCKTALK_ENABLE_ASAN OR CLOCKTALK_ENABLE_UBSAN)
  if(NOT (CMAKE_C_COMPILER_ID MATCHES "Clang" OR CMAKE_C_COMPILER_ID STREQUAL "GNU"))
    message(FATAL_ERROR "Sanitizers do not work reliably for compilers other than Clang or GCC.")
  endif()

  set(_SAN_FLAGS "")
  if(CLOCKTALK_ENABLE_ASAN)
    list(APPEND _SAN_FLAGS -fsanitize=address)
  endif()
  if(CLOCKTALK_ENABLE_UBSAN)
    list(APPEND _SAN_FLAGS -fsanitize=undefined)
    list(APPEND _SAN_FLAGS -fno-sanitize-recover=undefined)
  endif()
  list(APPEND _SAN_FLAGS -fno-omit-frame-pointer)

  add_compile_options(${_SAN_FLAGS})
  add_link_options(${_SAN_FLAGS})
endif()

add_compile_definitions(
    _POSIX_C_SOURCE=200809        # POSIX APIs (e.g., strdup, clock_gettime)
    CLOCKTALK_VERSION_MAJOR=1
    CLOCKTALK_VERSION_MINOR=0
    CLOCKTALK_VERSION_PATCH=0
    CLOCKTALK_NORELEASE           # indicate non-release build
)

include_directories(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)
