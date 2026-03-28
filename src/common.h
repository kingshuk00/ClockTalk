/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */

#ifndef CLOCKTALK_COMMON_H__
#define CLOCKTALK_COMMON_H__

#include<stdbool.h>

typedef struct {
  int diag;
  int error;
  bool timings;
  bool profile;
  bool pretty;
  bool opts;
} ClockTalkShowOpts;

typedef struct {
  struct {
    double len;
    int num;
    bool isOn;
  } win;

  struct {
    int rank;
    int num;
    bool isOn;
  } evt;
} ClockTalkMonOpts;

typedef struct {
  double eagerLimit;
  struct {
    bool overhead;
    bool flush;
    bool untraced;
  } ignore;
} ClockTalkSimOpts;

typedef struct ClockTalkOpts_type_ {
  char *filename;

  ClockTalkShowOpts show;

  ClockTalkMonOpts mon;

  ClockTalkSimOpts sim;

} ClockTalkOpts;

#endif  /* CLOCKTALK_COMMON_H__ */
