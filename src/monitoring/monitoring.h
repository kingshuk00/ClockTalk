/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */

#ifndef REPLAY_MONITORING_H__
#define REPLAY_MONITORING_H__

typedef struct ClockTalkOpts_type_ ClockTalkOpts;

extern void DoMonitoringEventBased(const ClockTalkOpts *const);
extern void DoMonitoringWindowed(const ClockTalkOpts *const);

#endif  /* REPLAY_MONITORING_H__ */
