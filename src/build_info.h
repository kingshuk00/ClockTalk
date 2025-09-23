/*
 * Copyright (c) 2025      High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef BUILD_INFO_H__
#define BUILD_INFO_H__

#define TOSTRING(s) STR(s)
#define STR(s) #s
#define VERSION_STR(m, n, p) TOSTRING(m)"."TOSTRING(n)"."TOSTRING(p)

#define CLOCKTALK_VERSION VERSION_STR(CLOCKTALK_VERSION_MAJOR, CLOCKTALK_VERSION_MINOR, CLOCKTALK_VERSION_PATCH)

#define CLOCKTALK_BUILD_DATE __DATE__
#define CLOCKTALK_BUILD_TIME __TIME__

#if defined(__clang__)
# define COMPILER_VERSION "Clang-"VERSION_STR(__clang_major__, __clang_minor__, __clang_patchlevel__)
#elif defined(__GNUC__)
# define COMPILER_VERSION "GCC-"VERSION_STR(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)
#else
# define COMPILER_VERSION "Unknown-x.y.z"
#endif

#endif  /* BUILD_INFO_H__ */
