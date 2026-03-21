/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Copyright (c) 2025      High Performance Computing Center Stuttgart,
 *                         University of Stuttgart. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */

#include"build_info.h"
#include"opt_parser.h"
#include"common.h"
#include"utils/utils.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#if !defined(__GLIBC__)
# error "Required `libc` is missing"
#endif
#include"argp.h"

const char *argp_program_version= "ClockTalk "CT_VERSION"\n  Compiled with "
  CT_COMPILER" on "CT_BUILD_DATE" at "CT_BUILD_TIME"\n\n"
  "  Copyright (C) 2026 Kingshuk Haldar. All rights reserved.\n\n"
  "  Copyright (C) 2025 High Performance Computing Center Stuttgart,\n"
  "                     University of Stuttgart. All rights reserved.\n";

const char *argp_program_bug_address=
  "<https://github.com/kingshuk00/ClockTalk/issues/new>";

inline static void interpretSpecialEvtsOpts(ClockTalkSimOpts *const opts,
                                            char *const optArg)
{
  char *ptr= strtok(optArg, ",\n ");
  while(NULL!= ptr) {
    if(0== strcmp("overhead", ptr)) {
      opts->ignore.overhead= true;
    } else if(0== strcmp("flush", ptr)) {
      opts->ignore.flush= true;
    } else if(0== strcmp("untraced", ptr)) {
      opts->ignore.untraced= true;
    } else {
      printf("Unknown process option for special events (%s)\n", ptr);
    }
    ptr= strtok(NULL, ",\n ");
  }
}

inline static void interpretEagerLimitOpt(ClockTalkSimOpts *const opts,
                                          char *const optArg)
{
  if(NULL!= optArg) {
    long x= atol(optArg);
    if(0== x) {
      Error("Invalid eager-limit specification \"%s\"\n", optArg);
      x= 32;
    }
    char u= 'k';
    sscanf(optArg, "%*d%c", &u);
    switch(u) {
    case 'B':
      opts->eagerLimit= (double) x;
      break;
    case 'M':
      opts->eagerLimit= (double)(x<<20);
      break;
    case 'G':
      opts->eagerLimit= (double)(x<<30);
      break;
    case 'k':                 /* fall-through */
    default:
      opts->eagerLimit= (double)(x<<10);
      break;
    }
  }
}

inline static void interpretMonTypes(ClockTalkMonOpts *const opts, char *const optArg)
{
  char *ptr= strtok(optArg, ",\n ");
  while(NULL!= ptr) {
    if(0== strcmp("window", ptr)) {
      opts->win.isOn= true;
    } else if(0== strcmp("event", ptr)) {
      opts->evt.isOn= true;
    } else {
      printf("Unknown process option for monitor-types (%s)\n", ptr);
    }
    ptr= strtok(NULL, ",\n ");
  }
}

static struct argp_option showOpts[]= {
  { "show-reviews", 'R', "1", OPTION_ARG_OPTIONAL, "Stepwise review level in stdout (default: 0)", 0 },
  { "show-errors", 'E', "1", OPTION_ARG_OPTIONAL, "Error level in stdout (default: 1)", 0 },
  { "show-timings", 'T', 0, 0, "I/O progress and timings in stdout (default: no)", 2 },
  { "export-profile", 'X', 0, 0, "Quick profile in a separate file (default: no)", 1 },
  { "pretty-output", 'P', 0, 0, "Formatted end-output in stdout (default: no)", 1 },
  { 0 }
};
static error_t parseShowOpts(int key, char *arg, struct argp_state *state)
{
  ClockTalkShowOpts *opts= state->input;
  switch(key) {
  case 'R':
    opts->diag= NULL!= arg? atoi(arg): 1;
    break;
  case 'E':
    opts->error= NULL!= arg? atoi(arg): 1;
    break;
  case 'T':
    opts->timings= true;
    break;
  case 'X':
    opts->profile= true;
    break;
  case 'P':
    opts->pretty= true;
    break;
  default:
    /* printf("Where are you (show)? (0x%x)\n", key); */
    return ARGP_ERR_UNKNOWN;
    break;
  }

  return 0;
}
static struct argp showOptsParser= { showOpts, parseShowOpts, 0 };

struct argp_option monOpts[]= {
  { "monitors", 'm', "window,event", 0, "Type of monitoring to perform", 0 },
  { "wmon-len", 2101, "1.0e9", 0, "Monitoring window in ns (default: 1e9 ns)", 1 },
  { "wmon-nwins", 2102, "1", 0, "#windows for simple moving average (default: 1)", 1 },
  { "emon-rank", 2201, "0", 0, "Event-based monitoring rank (default: 0)", 2 },
  { "emon-nevts", 2202, "1", 0, "#events accumulated per data-point (default: 1)", 2 },
  { 0 }
};
static error_t parseMonOpts(int key, char *arg, struct argp_state *state)
{
  ClockTalkMonOpts *opts= state->input;
  switch(key) {
  case 'm':
    interpretMonTypes(opts, arg);
    break;
  case 2101:
    opts->win.len= atof(arg);
    ErrorIf(opts->win.len< 0.9,
            "Invalid monitoring window length (%.9e ns)\n", opts->win.len);
    break;
  case 2102:
    opts->win.num= atoi(arg);
    ErrorIf(opts->win.num< 0.9,
            "Invalid #windows for moving-average (%d)\n", opts->win.num);
    break;
  case 2201:
    opts->evt.rank= atoi(arg);
    ErrorIf(opts->evt.rank< 0,
            "Invalid event-based monitoring rank (%d)\n", opts->evt.rank);
    break;
  case 2202:
    opts->evt.num= atoi(arg);
    ErrorIf(opts->evt.num< 0,
            "Invalid #events for event-based monitoring (%d)\n",
            opts->evt.num);
    break;
  default:
    /* printf("Where are you (mon)? (0x%x)\n", key); */
    return ARGP_ERR_UNKNOWN;
    break;
  }

  return 0;
}
static struct argp monOptsParser= { monOpts, parseMonOpts, 0 };

static struct argp_option simOpts[]= {
  { "eager-limit", 3001, "32k", 0, "Eager limit (default: 32k)" },
  { "ignore", 3002, "untraced,flush,overhead", 0, "Treat regions as useful (default: none)" },
  { 0 }
};
static error_t parseSimOpts(int key, char *arg, struct argp_state *state)
{
  ClockTalkSimOpts *opts= state->input;
  switch(key) {
  case 3001:
    interpretEagerLimitOpt(opts, arg);
    break;
  case 2102:
    interpretSpecialEvtsOpts(opts, arg);
    break;
  default:
    /* printf("Where are you (sim)? (0x%x)\n", key); */
    return ARGP_ERR_UNKNOWN;
    break;
  }

  return 0;
}
struct argp simOptsParser= { simOpts, parseSimOpts, 0 };

static struct argp_child childrenOpts[]= {
  { &showOptsParser, 0, "Display options:", 0 },
  { &monOptsParser, 0, "Mointor options:", 0 },
  { &simOptsParser, 0, "Calculation options:", 0 },
  { 0 }
};

static struct argp_option mainOpts[]= { { 0 } };
static error_t parseMainOpts(int, char *, struct argp_state *);

const char *const mainArgDesc= "<paraver-file-name>";
const char *const progDesc=
  "ClockTalk - Trace replay for critical path from Paraver trace files\n";
struct argp mainOptsParser= { mainOpts, parseMainOpts, mainArgDesc, progDesc, childrenOpts };
static error_t parseMainOpts(int key, char *arg, struct argp_state *state)
{
  ClockTalkOpts *const opts= (ClockTalkOpts *) state->input;
  switch(key) {
  case ARGP_KEY_ARG:
    FREE_IF(opts->filename);
    opts->filename= strdup(arg); /* strdup() is not C, but POSIX  */
    break;
  case ARGP_KEY_NO_ARGS:
    argp_usage(state);
    break;
  case ARGP_KEY_INIT:
    state->child_inputs[0]= &(opts->show);
    state->child_inputs[1]= &(opts->mon);
    state->child_inputs[2]= &(opts->sim);
    break;
  default:
    /* printf("Where are you (main)? (0x%x)\n", key); */
    return ARGP_ERR_UNKNOWN;
    break;
  }

  return 0;
}

ClockTalkOpts *ParseOpts(const int argc, char **argv)
{
  ClockTalkOpts *opts= (ClockTalkOpts *) malloc(sizeof(ClockTalkOpts));
  memset(opts, 0, sizeof(ClockTalkOpts));

  argp_parse(&mainOptsParser, argc, argv, 0, 0, opts);

  if(NULL== opts->filename) {
    goto bad;
  }

  if(opts->mon.win.isOn) {
    if(opts->mon.win.len< 0.9) {
      printf("Windowed monitoring: window-length is invalid (1.0e9 ns)\n");
      opts->mon.win.len= 1.0e9;
    }
#if 0
    if(opts->mon.win.num< 1) {
      printf("Windowed monitoring: #windows for SMA is invalid (1)\n");
      opts->mon.win.num= 1;
    }
#endif
  }

  if(opts->mon.evt.isOn) {
    if(opts->mon.evt.rank< 0) {
      printf("Event-based monitoring: rank is invalid (0)\n");
      opts->mon.evt.rank= 0;
    }
    if(opts->mon.evt.num< 1) {
      printf("Event-based monitoring: #events per report is invalid (1)\n");
      opts->mon.evt.num= 1;
    }
  }

  UtilSetShowFunctions(opts->show.error, opts->show.diag);

  goto bye;

 bad:
  FREE_IF(opts);
  argp_help(&mainOptsParser, stdout, ARGP_HELP_LONG, NULL);

 bye:
  return opts;
}

