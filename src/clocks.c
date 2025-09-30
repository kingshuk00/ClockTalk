/*
 * Copyright (c) 2024-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#include"common.h"
#include"utils.h"
#include"paraver.h"
#include"clocks.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

inline static bool isPaused(const int s) { return 0!= s; }
inline static bool isPlaying(const int s) { return 0== s; }
inline static bool isDisabled(const int s) { return -2== s; }

inline static const char *name(const int e) { return GetParaverMPIEvtName(e); }

static struct {
  int *state;
  double *since;
  double *onSince;   /* positive: tracing on; negative: tracing off */
} curr;

inline static void setState(const int p, const int s) { curr.state[p]= s; }
inline static void setSince(const int p, const double t) { curr.since[p]= t; }
inline static void setCurr(const int p, const double at, const int id) { setState(p, id); setSince(p, at); }

inline static int state(const int p) { return curr.state[p]; }
inline static bool isCurrPaused(const int p) { return isPaused(state(p)); }
inline static bool isCurrPlaying(const int p) { return isPlaying(state(p)); }
inline static bool isCurrDisabled(const int p) { return isDisabled(state(p)); }

inline static const char *currName(const int p) { return name(state(p)); }

inline static double since(const int p) { return curr.since[p]; }

inline static void setTracing(const int p, const double t) { curr.onSince[p]= t; }
inline static double tracingSince(const int p) { return curr.onSince[p]; }
inline static bool isTracing(const int p) { return tracingSince(p)> -0.1; }

inline static double enabledAt(const int p) { return isTracing(p)? tracingSince(p): -1.0; }
inline static double disabledAt(const int p) { return isTracing(p)? -1.0: -tracingSince(p); }
inline static void disable(const int p, const double at)
{
  ErrorIf(!isTracing(p),
          "%d: disabled since %.0lf, overwritten at %.0lf\n", p, disabledAt(p), at);
  setTracing(p, -at);
}
inline static void enable(const int p, const double at)
{
  ErrorIf(isTracing(p)&& !SameTime(enabledAt(p), at),
          "%d: enabled since %.0lf, overwritten at %.0lf\n", p, enabledAt(p), at);
  setTracing(p, at);
}

ClockType Clocks;

inline static double elapsed(const int p) { return Clocks.elapsed[p]; }
inline static void setElapsed(const int p, const double t) { Clocks.elapsed[p]= t; }

inline static void setTraced(const int p, const double t) { Clocks.traced[p]= t; }
inline static void updateTraced(const int p, const double t) { Clocks.traced[p]+= t- since(p); }

inline static void updateFlush(const int p, const double t) { Clocks.flush[p]+= t- since(p); }

inline static void addUseful(const int p, const double t) { Clocks.useful[p]+= t; }
inline static void updateUseful(const int p, const double t) { Clocks.useful[p]+= t- since(p); }

inline static double critical(const int p) { return Clocks.critical[p]; }
inline static void setCritical(const int p, const double t) { Clocks.critical[p]= t; }
inline static void updateCritical(const int p, const double t) { Clocks.critical[p]+= t- since(p); }

inline static void elapse(const int p, const double t, const int evt) { setCurr(p, t, evt); setElapsed(p, t); }

inline static void currCalc(const int p, const double t, const int evt)
{
  switch(state(p)) {
  case 0:                       /* useful */
    if(isTracing(p)) {
      updateCritical(p, t);
      updateUseful(p, t);
    }
    break;
  case -2:                      /* disabled */
    Debug1("%d: %s at %.0lf in disabled state\n", p, name(evt), t);
    break;
  case -3:                      /* flush */
    updateFlush(p, t);
    break;
  case -99:                     /* invalid */
  case -4:                      /* trace-init */
    break;
  case -1:                      /* ended */
    Error("%d: %s -> %s at %.0lf\n", p, currName(p), name(evt), t);
    break;
  default:                      /* MPI */
    break;
  }

  if(isTracing(p)) {
    updateTraced(p, t);
  }
}

/* evt should always be zero here - starting a useful region */
void ClockPlay(const int p, const double t, const int evt)
{
  ErrorIf(0!= evt, "%d: play clock at %.0lf with event %d\n", p, t, evt);

  currCalc(p, t, evt);

  if(!isTracing(p)) {
    enable(p, t);
  }

  elapse(p, t, evt);

  if(GlOpts.evt_mon.rank== p) {
    /* /\* TlOutput(maxElapsed(), maxTraced(p), Clocks.critical[pOfMaxCritical()], maxUseful(), avgUseful()); *\/ */
    /* TlOutput(elapsed(p), traced(p), critical(p), useful(p), avgUseful()); */
  }

  Debug1("%d: clock-play at %.0lf, critical: %.0lf\n", p, elapsed(p),
         critical(p));
}

/* evt should always be MPI (> 0) here - starting an MPI region */
void ClockPauseMPI(const int p, const double t, const int evt)
{
  ErrorIf(!isTracing(p),
          "%d: pause clock event (%d) at %.0lf, tracing disabled at %.0lf\n", p,
          evt, t, disabledAt(p));
  ErrorIf(isCurrPaused(p)&& !SameTime(t, since(p)),
          "%d: pause clock event (%d) at %.0lf, paused since %.0lf\n", p,
          evt, t, since(p));

  currCalc(p, t, evt);

  elapse(p, t, evt);

  Debug1("%d: clock-pause(%s) at %.0lf, critical %.0lf\n", p, currName(p),
         elapsed(p), critical(p));
}
/* evt shoudl be < 0 here */
void ClockPauseTrace(const int p, const double t, const int evt)
{
  currCalc(p, t, evt);

  switch(evt) {
  case -99:
    Error("%d: %s -> %s at %.0lf\n", p, currName(p), name(evt), t);
    break;
  case -1:
    Debug1("%d: end tracing at %.0lf\n", p, t);
    break;
  case -2:
    if(isTracing(p)) {
      disable(p, t);
    }
    break;
  case -3:
    Debug1("%d: start flush at %.0lf\n", p, t);
    break;
  case -4:
    Debug1("%d: start trace-init at %.0lf\n", p, t);
    break;
  default:
    Error("%d: %s -> unknown(%d) at %.0lf\n", p, currName(p), evt, t);
    break;
  }

  elapse(p, t, evt);

  Debug1("%d: clock-pause(%s) at %.0lf, critical %.0lf\n", p, currName(p),
         elapsed(p), critical(p));
}

void ClockStart(const int np, const double t0, const double *const pt0s)
{
  for(int ip= 0; ip< np; ++ip) {
    const int pt0= pt0s[ip* 2];

    enable(ip, t0);
    setCritical(ip, pt0);
    setTraced(ip, pt0);
    elapse(ip, pt0, 0);
  }

  /* TlCreateFile(np); */
}
void ClockEnd(const int p, const double t, const double universeEnds)
{
  /* even if tracing is disabled, this is performed */
  if(isCurrPlaying(p)) {
    updateUseful(p, t);
  } else if(isCurrDisabled(p)) {
    addUseful(p, universeEnds- t);
  }
  updateCritical(p, t);
  updateTraced(p, t);
  elapse(p, t, -1);
  /* /\* TlOutput(elapsed(p), traced(p), critical(p), useful(p)); *\/ */
  /* TlOutput(maxElapsed(), maxTraced(), Clocks.critical[pOfMaxCritical()], maxUseful(), avgUseful()); */
  /* TlCloseFile(); */
  Debug1("%d: end clock at %.0lf, critical %.0lf\n", p, elapsed(p), critical(p));
}

void ClockInit(const int np)
{
  Clocks.elapsed= (double *) malloc(sizeof(double)* np);
  memset(Clocks.elapsed, 0, sizeof(double)* np);

  Clocks.traced= (double *) malloc(sizeof(double)* np);
  memset(Clocks.traced, 0, sizeof(double)* np);

  Clocks.flush= (double *) malloc(sizeof(double)* np);
  memset(Clocks.flush, 0, sizeof(double)* np);

  Clocks.useful= (double *) malloc(sizeof(double)* np);
  memset(Clocks.useful, 0, sizeof(double)* np);

  Clocks.critical= (double *) malloc(sizeof(double)* np);
  memset(Clocks.critical, 0, sizeof(double)* np);

  curr.state= (int *) malloc(sizeof(int)* np);
  memset(curr.state, 0, sizeof(int)* np);

  curr.since= (double *) malloc(sizeof(double)* np);
  memset(curr.since, 0, sizeof(double)* np);

  curr.onSince= (double *) malloc(sizeof(double)* np);
  memset(curr.onSince, 0, sizeof(double)* np);

  Clocks.np= np;
}
void ClockFinalize()
{
  FREE_IF(curr.onSince);
  FREE_IF(curr.since);
  FREE_IF(curr.state);

  FREE_IF(Clocks.critical);
  FREE_IF(Clocks.useful);
  FREE_IF(Clocks.flush);
  FREE_IF(Clocks.traced);
  FREE_IF(Clocks.elapsed);
}
