/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */

#ifndef CLOCKTALK_REPLAY_REPLAY_H__
#define CLOCKTALK_REPLAY_REPLAY_H__


#include<stdbool.h>

typedef struct ClockTalkOpts_type_ ClockTalkOpts;

extern void ReplayTrace(const ClockTalkOpts *const);
extern void ShowStatsPostReplay(const bool);


#endif  /* CLOCKTALK_REPLAY_REPLAY_H__ */
