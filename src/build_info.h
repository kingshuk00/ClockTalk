/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Copyright (c) 2025      High Performance Computing Center Stuttgart,
 *                         University of Stuttgart. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */

#ifndef CLOCKTALK_BUILD_INFO_H__
#define CLOCKTALK_BUILD_INFO_H__

#define TOSTRING(s) STR(s)
#define STR(s) #s
#define VERSION_STR(m, n, p) TOSTRING(m)"."TOSTRING(n)"."TOSTRING(p)

#if defined(CLOCKTALK_NORELEASE)
# define CT_VERSION_DEVEL x
#else
# define CT_VERSION_DEVEL
#endif
#define CT_VERSION VERSION_STR(CLOCKTALK_VERSION_MAJOR, CLOCKTALK_VERSION_MINOR, CLOCKTALK_VERSION_PATCH)TOSTRING(CT_VERSION_DEVEL)

#define CT_BUILD_DATE __DATE__
#define CT_BUILD_TIME __TIME__

#if defined(__INTEL_LLVM_COMPILER)
# define CT_COMPILER __VERSION__
#elif defined(__clang__)
# define CT_COMPILER "Clang-"VERSION_STR(__clang_major__, __clang_minor__, __clang_patchlevel__)
#elif defined(__GNUC__)
# define CT_COMPILER "GCC-"VERSION_STR(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)
#else
# define CT_COMPILER "Unknown"
#endif

#endif  /* CLOCKTALK_BUILD_INFO_H__ */
