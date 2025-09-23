/*
 * Copyright (c) 2024-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef REPLAY_CLOCKS_H__
#define REPLAY_CLOCKS_H__


extern void ClockInit(const int);
extern void ClockFinalize();

extern void ClockPlay(const int, const double, const int);
extern void ClockPauseMPI(const int, const double, const int);
extern void ClockPauseTrace(const int, const double, const int);

extern void ClockStart(const int, const double, const double *const);
extern void ClockEnd(const int, const double, const double);

typedef struct {
  double *elapsed;              /* also runs when disabled */
  double *traced;               /* runs unless disabled */
  double *flush;                /* runs when flushing */
  double *useful;               /* runs for useful events */
  double *critical;             /* runs for useful and wait */
  int np;
} ClockType;

extern ClockType Clocks;

/* state-identifiers:
 *  0: computation aka useful
 * >0: MPI
 * -1: tracing ended
 * -2: tracing disabled - ideal includes, but not useful
 * -3: flush
 * -4: trace-init
 * -99: invalid
 */


inline static double ClockGetElapsed(const int p) { return Clocks.elapsed[p]; }
inline static double ClockGetTraced(const int p) { return Clocks.traced[p]; }
inline static double ClockGetUseful(const int p) { return Clocks.useful[p]; }
inline static double ClockGetCritical(const int p) { return Clocks.critical[p]; }

inline static void ClockDebug(const int p)
{
  printf("Clocks(%d): elapsed= %.0lf, traced= %.0lf, useful= %.0lf, critical= %.0lf\n",
         p,
         ClockGetElapsed(p), ClockGetTraced(p), ClockGetUseful(p), ClockGetCritical(p));
}

inline static void ClockSetCritical(const int p, const double t)
{
  Debug1("%d: set critical at %.0lf: %.0lf -> %.0lf\n", p,
         ClockGetElapsed(p), ClockGetCritical(p), t);
  Clocks.critical[p]= t;
}
inline static void ClockUpdateCritical(const int p, const double delta_t)
{
  Debug1("%d: update critical at %.0lf: %.0lf -> %.0lf\n", p,
         ClockGetElapsed(p), ClockGetCritical(p), ClockGetCritical(p)+ delta_t);
  Clocks.critical[p]+= delta_t;
}

inline static double ClockGetMaxElapsed(const int np)
{
  double elapsed= 0.0;
  for(int ip= 0; ip< np; ++ip) {
    elapsed= MAX(ClockGetElapsed(ip), elapsed);
  }
  return elapsed;
}

inline static double ClockGetMaxTraced(const int np)
{
  double traced= 0.0;
  for(int ip= 0; ip< np; ++ip) {
    traced= MAX(ClockGetTraced(ip), traced);
  }
  return traced;
}

inline static double ClockGetMaxCritical(const int np)
{
  double critical= 0.0;
  for(int ip= 0; ip< np; ++ip) {
    critical= MAX(ClockGetCritical(ip), critical);
  }
  return critical;
}

inline static double ClockGetMaxUseful(const int np)
{
  double useful= 0.0;
  for(int ip= 0; ip< np; ++ip) {
    useful= MAX(ClockGetUseful(ip), useful);
  }
  return useful;
}

inline static double ClockGetAvgUseful(const int np)
{
  double avg= 0.0;
  for(int ip= 0; ip< np; ++ip) {
    avg+= ClockGetUseful(ip);
  }
  return avg/ ((double) np);
}

#endif  /* REPLAY_CLOCKS_H__  */
