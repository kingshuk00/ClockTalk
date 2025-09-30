/*
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef REPLAY_COMMON_H__
#define REPLAY_COMMON_H__
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>
#include<time.h>

typedef struct {
  char *filename;

  struct {
    int diag;
    int error;
    bool timings;
    bool profile;
    bool pretty;
  } show_opts;

  struct {
    double win_len;
    int nwins_sma;
    bool enabled;
  } win_mon;

  struct {
    int rank;
    int nevts_report;
    bool enabled;
  } evt_mon;

  struct {
    double eager_limit;
    struct {
      bool trace_evts;
      bool flush_evts;
      bool disabled_tracing;
    } ignore;
  } sim_opts;
} GlobalOpts;

extern GlobalOpts GlOpts;


#endif  /* REPLAY_COMMON_H__ */
