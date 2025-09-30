/*
 * Copyright (c) 2025      High Performance Computing Center Stuttgart,
 *                         University of Stuttgart. All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#include"build_info.h"
#include"arg_opt_parser.h"
#include"common.h"
#include"utils.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#if !defined(__GLIBC__)
# error "Required `libc` is missing"
#endif
#include"argp.h"

static void printVersion(const char *const name)
{
  printf("%s (%s on %s at %s) %s\n\n", name, COMPILER_VERSION,
         CLOCKTALK_BUILD_DATE, CLOCKTALK_BUILD_TIME, CLOCKTALK_VERSION);
  printf("Copyright (C) 2025 High Performance Computing Center Stuttgart,\n"
         "                   University of Stuttgart. All rights reserved.\n");
}

GlobalOpts GlOpts= { NULL, { 0, 1, false, false, false }, { 0.0, -1, false }, { -1, 0, false }, {32768.0, { false, false, false } } };

inline static void interpretSpecialEvtsOpts(GlobalOpts *const opts,
                                            char *const optArg)
{
  char *ptr= strtok(optArg, ",\n ");
  while(NULL!= ptr) {
    if(0== strcmp("overhead", ptr)) {
      opts->sim_opts.ignore.trace_evts= true;
    } else if(0== strcmp("flush", ptr)) {
      opts->sim_opts.ignore.flush_evts= true;
    } else if(0== strcmp("traceability", ptr)) {
      opts->sim_opts.ignore.disabled_tracing= true;
    } else {
      printf("Unknown process option for special events (%s)\n", ptr);
    }
    ptr= strtok(NULL, ",\n ");
  }
}

inline static void interpretEagerLimitOpt(GlobalOpts *const opts,
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
      opts->sim_opts.eager_limit= (double) x;
      break;
    case 'M':
      opts->sim_opts.eager_limit= (double)(x<<20);
      break;
    case 'G':
      opts->sim_opts.eager_limit= (double)(x<<30);
      break;
    case 'k':                 /* fall-through */
    default:
      opts->sim_opts.eager_limit= (double)(x<<10);
      break;
    }
  } else {
    /* printf("Problem\n"); */
  }
}

inline static void interpretMonTypes(GlobalOpts *const opts, char *const optArg)
{
  char *ptr= strtok(optArg, ",\n ");
  while(NULL!= ptr) {
    if(0== strcmp("window", ptr)) {
      opts->win_mon.enabled= true;
    } else if(0== strcmp("event", ptr)) {
      opts->evt_mon.enabled= true;
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
  GlobalOpts *opts= state->input;
  switch(key) {
  case 'R':
    opts->show_opts.diag= NULL!= arg? atoi(arg): 1;
    break;
  case 'E':
    opts->show_opts.error= NULL!= arg? atoi(arg): 1;
    break;
  case 'T':
    opts->show_opts.timings= true;
    break;
  case 'X':
    opts->show_opts.profile= true;
    break;
  case 'P':
    opts->show_opts.pretty= true;
    break;
  case ARGP_KEY_ARG:
    printf("ARGP_KEY_ARG(show)\n");
    break;
  case ARGP_KEY_ARGS:
    printf("ARGP_KEY_ARGS(show)\n");
    break;
  case ARGP_KEY_NO_ARGS:
    /* printf("ARGP_KEY_NO_ARGS(show)\n"); */
    break;
  case ARGP_KEY_INIT:
    /* printf("ARGP_KEY_INIT(show)\n"); */
    break;
  case ARGP_KEY_END:
    /* printf("ARGP_KEY_END(show)\n"); */
    break;
  case ARGP_KEY_SUCCESS:
    /* printf("ARGP_KEY_SUCCESS(show)\n"); */
    break;
  case ARGP_KEY_FINI:
    /* printf("ARGP_KEY_FINI(show)\n"); */
    break;
  case ARGP_KEY_ERROR:
    printf("ARGP_KEY_ERROR(show)\n");
    break;
  default:
    printf("Where are you (show)? (0x%x)\n", key);
    return ARGP_ERR_UNKNOWN;
    break;
  }

  return 0;
}
static struct argp showOptsParser= { showOpts, parseShowOpts, 0 };

struct argp_option monOpts[]= {
  { "monitors", 'm', "window,event", 0, "Type of monitoring to perform", 0 },
  { "wmon-len", 2101, "1.0e9", 0, "Monitoring window in ns (default: 1e9 ns)", 1 },
  { "wmon-sma", 2102, "1", 0, "#windows for simple moving average (default: 1)", 1 },
  { "emon-rank", 2201, "0", 0, "Event-based monitoring rank (default: 0)", 2 },
  { "emon-nevts", 2202, "1", 0, "#events accumulated per data-point (default: 1)", 2 },
  { 0 }
};
static error_t parseMonOpts(int key, char *arg, struct argp_state *state)
{
  GlobalOpts *opts= state->input;
  switch(key) {
  case 'm':
    interpretMonTypes(opts, arg);
    break;
  case 2101:
    opts->win_mon.win_len= atof(arg);
    ErrorIf(opts->win_mon.win_len< 0.9,
            "Invalid monitoring window length (%.9e ns)\n", opts->win_mon.win_len);
    break;
  case 2102:
    opts->win_mon.nwins_sma= atoi(arg);
    ErrorIf(opts->win_mon.nwins_sma< 0.9,
            "Invalid #windows for moving-average (%d)\n", opts->win_mon.nwins_sma);
    break;
  case 2201:
    opts->evt_mon.rank= atoi(arg);
    ErrorIf(opts->evt_mon.rank< 0,
            "Invalid event-based monitoring rank (%d)\n", opts->evt_mon.rank);
    break;
  case 2202:
    opts->evt_mon.nevts_report= atoi(arg);
    ErrorIf(opts->evt_mon.nevts_report< 0,
            "Invalid #events for event-based monitoring (%d)\n",
            opts->evt_mon.nevts_report);
    break;
  case ARGP_KEY_ARG:
    printf("ARGP_KEY_ARG(mon)\n");
    break;
  case ARGP_KEY_ARGS:
    printf("ARGP_KEY_ARGS(mon)\n");
    break;
  case ARGP_KEY_NO_ARGS:
    /* printf("ARGP_KEY_NO_ARGS(mon)\n"); */
    break;
  case ARGP_KEY_INIT:
    /* printf("ARGP_KEY_INIT(mon)\n"); */
    break;
  case ARGP_KEY_END:
    /* printf("ARGP_KEY_END(mon)\n"); */
    break;
  case ARGP_KEY_SUCCESS:
    /* printf("ARGP_KEY_SUCCESS(mon)\n"); */
    break;
  case ARGP_KEY_FINI:
    /* printf("ARGP_KEY_FINI(mon)\n"); */
    break;
  case ARGP_KEY_ERROR:
    printf("ARGP_KEY_ERROR(mon)\n");
    break;
  default:
    printf("Where are you (mon)? (0x%x)\n", key);
    return ARGP_ERR_UNKNOWN;
    break;
  }

  return 0;
}
static struct argp monOptsParser= { monOpts, parseMonOpts, 0 };

static struct argp_option simOpts[]= {
  { "eager-limit", 3001, "32k", 0, "Eager limit (default: 32k)" },
  { "ignore-events", 3002, "traceability,flush,overhead", 0, "Trace-events as useful (default: none)" },
  { 0 }
};
static error_t parseSimOpts(int key, char *arg, struct argp_state *state)
{
  GlobalOpts *opts= state->input;
  switch(key) {
  case 3001:
    interpretEagerLimitOpt(opts, arg);
    break;
  case 2102:
    interpretSpecialEvtsOpts(opts, arg);
    break;
  case ARGP_KEY_ARG:
    printf("ARGP_KEY_ARG(sim)\n");
    break;
  case ARGP_KEY_ARGS:
    printf("ARGP_KEY_ARGS(sim)\n");
    break;
  case ARGP_KEY_NO_ARGS:
    /* printf("ARGP_KEY_NO_ARGS(sim)\n"); */
    break;
  case ARGP_KEY_INIT:
    /* printf("ARGP_KEY_INIT(sim)\n"); */
    break;
  case ARGP_KEY_END:
    /* printf("ARGP_KEY_END(sim)\n"); */
    break;
  case ARGP_KEY_SUCCESS:
    /* printf("ARGP_KEY_SUCCESS(sim)\n"); */
    break;
  case ARGP_KEY_FINI:
    /* printf("ARGP_KEY_FINI(sim)\n"); */
    break;
  case ARGP_KEY_ERROR:
    printf("ARGP_KEY_ERROR(sim)\n");
    break;
  default:
    printf("Where are you (sim)? (0x%x)\n", key);
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

static struct argp_option mainOpts[]= {
  {"version", 'V', NULL, 0, NULL, 0 },
  { 0 }
};
static error_t parseMainOpts(int, char *, struct argp_state *);

const char *const mainArgDesc= "<paraver-file-name>";
const char *const progDesc=
  "ClockTalk - a trace replay tool for critical path from Paraver file\n\n"
  "Program options:";
struct argp mainOptsParser= { mainOpts, parseMainOpts, mainArgDesc, progDesc, childrenOpts };
static error_t parseMainOpts(int key, char *arg, struct argp_state *state)
{
  GlobalOpts *opts= state->input;
  switch(key) {
  case 'V':
    printVersion(state->name);
    exit(0);
    break;
  case ARGP_KEY_ARG:
    opts->filename= strdup(arg); /* this is not C, but POSIX  */
    break;
  case ARGP_KEY_ARGS:
    printf("ARGP_KEY_ARGS(main)\n");
    break;
  case ARGP_KEY_NO_ARGS:
    argp_usage(state);
    break;
  case ARGP_KEY_INIT:
    for(int i= 0; i< 3; ++i) {
      state->child_inputs[i]= &GlOpts;
    }
    /* printf("ARGP_KEY_INIT(main)\n"); */
    break;
  case ARGP_KEY_END:
    /* printf("ARGP_KEY_END(main)\n"); */
    break;
  case ARGP_KEY_SUCCESS:
    /* printf("ARGP_KEY_SUCCESS(main)\n"); */
    break;
  case ARGP_KEY_FINI:
    /* printf("ARGP_KEY_FINI(main)\n"); */
    break;
  case ARGP_KEY_ERROR:
    printf("ARGP_KEY_ERROR(main)\n");
    break;
  default:
    printf("Where are you (main)? (0x%x)\n", key);
    return ARGP_ERR_UNKNOWN;
    break;
  }

  return 0;
}

int ParseArgs(const int argc, char **argv)
{
  argp_parse(&mainOptsParser, argc, argv, 0, 0, &GlOpts);
  int ret= 0;

  if(NULL== GlOpts.filename) {
    ret= 1;
  }

  if(GlOpts.win_mon.enabled) {
    if(GlOpts.win_mon.win_len< 0.9) {
      printf("Windowed monitoring: window-length is invalid (1.0e9 ns)\n");
      GlOpts.win_mon.win_len= 1.0e9;
    }
#if 0
    if(GlOpts.win_mon.nwins_sma< 1) {
      printf("Windowed monitoring: #windows for SMA is invalid (1)\n");
      GlOpts.win_mon.nwins_sma= 1;
    }
#endif
  }

  if(GlOpts.evt_mon.enabled) {
    if(GlOpts.evt_mon.rank< 0) {
      printf("Event-based monitoring: rank is invalid (0)\n");
      GlOpts.evt_mon.rank= 0;
    }
    if(GlOpts.evt_mon.nevts_report< 1) {
      printf("Event-based monitoring: #events per report is invalid (1)\n");
      GlOpts.evt_mon.nevts_report= 1;
    }
  }

  return ret;
}

void ArgHelp()
{
  argp_help(&mainOptsParser, stdout, ARGP_HELP_LONG, NULL);
}
