/*
 * Copyright (c) 2025      High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#include"monitoring.h"
#include"trace_data.h"
#include"common.h"
#include<float.h>

inline static FILE *emFileOpen()
{
  const int len= strlen(GlOpts.filename)+ 8;
  char *fn= (char *) malloc(sizeof(char)* len);
  memset(fn, 0, sizeof(char)* len);
  memcpy(fn, GlOpts.filename, sizeof(char)* (len- 8));
  /* printf("fn0= \"%s\"\n", fn); */
  strcat(fn, ".em.dat");
  /* printf("fn1= \"%s\", fn[-3]= '%c', fn[-2]= '%c', fn[-1]= '%c'\n", fn, fn[len- 3], fn[len- 2], fn[len- 1]); */

  FILE *fp= fopen(fn, "w");
  FREE_IF(fn);
  return fp;
}

#define EM_TRAC 0
#define EM_CRIT 1
#define EM_UMAX 2
#define EM_UAVG 3
#define EM_LAST EM_UAVG
#define EM_NCOMPS (EM_LAST+1)

static struct {
  FILE *fp;
  double *history[EM_NCOMPS];
  int nspans;

  struct {
    double elapsed;
    double *onSince;
    double traced;
    double critical;
    double *useful;
    double pfactor;
    struct {
      double *since;
      int *state;
    } curr;
    int np;
  } clocks;
} em= { NULL, { NULL }, 2, { 0.0, NULL, 0.0, 0.0, NULL, 0.0, { NULL, NULL }, 0 } };
inline static void emPushHistory(const double traced, const double critical,
                                 const double umax, const double uavg)
{
  for(int i= EM_TRAC; i< EM_NCOMPS; ++i) {
    memmove(em.history[i], em.history[i]+ 1, sizeof(double)* em.nspans);
  }
  em.history[EM_TRAC][em.nspans]= traced;
  em.history[EM_CRIT][em.nspans]= critical;
  em.history[EM_UMAX][em.nspans]= umax;
  em.history[EM_UAVG][em.nspans]= uavg;
}
inline static void emUsefuls(double *const umaxp, double *const uavgp)
{
  double umax= 0.0, uavg= 0.0;
  for(int ip= 0; ip< em.clocks.np; ++ip) {
    uavg+= em.clocks.useful[ip];
    umax= MAX(em.clocks.useful[ip], umax);
  }
  uavg*= em.clocks.pfactor;
  *umaxp= umax; *uavgp= uavg;
}
inline static double emTrfCum()
{
  return em.history[EM_TRAC][em.nspans]- em.history[EM_CRIT][em.nspans]> 1.0e-12?
         em.history[EM_CRIT][em.nspans]/ em.history[EM_TRAC][em.nspans]: 1.0;
}
inline static double emSerCum()
{
  return em.history[EM_CRIT][em.nspans]- em.history[EM_UMAX][em.nspans]> 1.0e-12?
         em.history[EM_UMAX][em.nspans]/ em.history[EM_CRIT][em.nspans]: 1.0;
}
inline static double emLbeCum()
{
  return em.history[EM_UMAX][em.nspans]- em.history[EM_UAVG][em.nspans]> 1.0e-12?
         em.history[EM_UAVG][em.nspans]/ em.history[EM_UMAX][em.nspans]: 1.0;
}
inline static double emTrfLoc()
{
  const double t= em.history[EM_TRAC][em.nspans]- em.history[EM_TRAC][0];
  const double c= em.history[EM_CRIT][em.nspans]- em.history[EM_CRIT][0];
  return t- c> 1.0e-12? c/ t: 1.0;
}
inline static double emSerLoc()
{
  const double c= em.history[EM_CRIT][em.nspans]- em.history[EM_CRIT][0];
  const double m= em.history[EM_UMAX][em.nspans]- em.history[EM_UMAX][0];
  return c- m> 1.0e-12? m/ c: 1.0;
}
inline static double emLbeLoc()
{
  const double m= em.history[EM_UMAX][em.nspans]- em.history[EM_UMAX][0];
  const double a= em.history[EM_UAVG][em.nspans]- em.history[EM_UAVG][0];
  /* if(m< a) { m= a; } */
  return m- a> 1.0e-12? a/ m: 1.0;
}
inline static void emOutput()
{
  /* input times are in ns */
  const double factor= 1.0e-9;  /* ns to s */
  const double elapsed= em.clocks.elapsed* factor;
  const double traced= em.clocks.traced* factor;
  const double critical= em.clocks.critical* factor;
  double umax= 0.0, uavg= 0.0; emUsefuls(&umax, &uavg);
  umax*= factor; uavg*= factor;

  fprintf(em.fp, "%.9e %.9e %.9e %.9e %.9e", elapsed, traced, critical, umax,
          uavg); /* 1 2 3 4 5 */

  emPushHistory(traced, critical, umax, uavg);
  fprintf(em.fp, " %.9e %.9e %.9e", emLbeCum(), emSerCum(),
          emTrfCum()); /* 6 7 8 */
  fprintf(em.fp, " %.9e %.9e %.9e\n", emLbeLoc(), emSerLoc(),
          emTrfLoc()); /* 9 10 11 */
}
inline static void emInit(const int np)
{
  em.fp= emFileOpen();
  fprintf(em.fp, "%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n",
          "#elapsed-1", "traced-2", "ideal-3",
          "max-useful-4", "avg-useful-5", "cum-load-bal-6", "cum-ser-eff-7",
          "cum-xfer-eff-8",
          "loc-load-bal-9", "loc-ser-eff-10", "loc-xfer-eff-11");
  em.nspans= GlOpts.evt_mon.nevts_report;
  em.history[0]= (double *) malloc(sizeof(double)* (em.nspans+ 1)* EM_NCOMPS);
  memset(em.history[0], 0, sizeof(double)* (em.nspans+ 1)* EM_NCOMPS);
  for(int i= 1; i< EM_NCOMPS; ++i) {
    em.history[i]= em.history[i- 1]+ (em.nspans+ 1);
  }

  em.clocks.elapsed= 0.0;
  em.clocks.traced= 0.0;
  em.clocks.critical= 0.0;

  em.clocks.np= np;
  em.clocks.pfactor= 1.0/ ((double) np);

  em.clocks.onSince= (double *) malloc(sizeof(double)* np);
  memset(em.clocks.onSince, 0, sizeof(double)* np);

  em.clocks.useful= (double *) malloc(sizeof(double)* np);
  memset(em.clocks.useful, 0, sizeof(double)* np);

  em.clocks.curr.since= (double *) malloc(sizeof(double)* np);
  memset(em.clocks.curr.since, 0, sizeof(double)* np);

  em.clocks.curr.state= (int *) malloc(sizeof(int)* np);
  memset(em.clocks.curr.state, 0, sizeof(int)* np);

  emOutput();
}
inline static void emFinalize()
{
  FREE_IF(em.clocks.curr.state);
  FREE_IF(em.clocks.curr.since);
  FREE_IF(em.clocks.useful);
  FREE_IF(em.clocks.onSince);
  FREE_IF(em.history[0]);
  memset(em.history, 0, sizeof(double *)* EM_NCOMPS);
  if(NULL!= em.fp) {
    fclose(em.fp); em.fp= NULL;
  }
}

#undef EM_NCOMPS
#undef EM_LAST
#undef EM_UAVG
#undef EM_UMAX
#undef EM_CRIT
#undef EM_TRAC

void DoMonitoringEventBased()
{
  const int np= TraceGetNumProcs();
  emInit(np);
  const int r= GlOpts.evt_mon.rank;
  double *const useful= em.clocks.useful;
  int *const state= em.clocks.curr.state;
  double *const since= em.clocks.curr.since;
  for(TraceResetIterEvts(); TraceRemainEvts(); TraceIncrIterEvts()) {
    const double t= TraceGetCurrEvtAt();
    const int e= TraceGetCurrEvtId();
    const int p= TraceGetCurrEvtProc();
    if(p== r) {                 /* useful + critical for this rank */
      const double delt= t- since[p];
      em.clocks.elapsed+= delt;
      if(em.clocks.onSince[p]> -0.1) {
        em.clocks.traced+= delt;
      }

      if(0== e) {               /* something (MPI/special) -> useful, this is a point */
        if(TraceGetCurrEvtCrit()> 0.1) {
          em.clocks.critical= TraceGetCurrEvtCrit();
        }
        for(int ip= 0; ip< np; ++ip) {
          if(0== state[ip]) {
            useful[ip]+= t- since[ip];
          }
          since[ip]= t;
        }
        emOutput();
      } else {                  /* useful -> something(MPI/special) */
        if(true|| 0== state[p]) {
          useful[p]+= delt;
          em.clocks.critical+= delt;
        }
        since[p]= t;
        state[p]= e;
      }
    } else {                    /* just usefuls for these ranks */
      if(0== state[p]) {
        useful[p]+= t- since[p];
      }
      since[p]= t;
      state[p]= e;
    }
  }
  emFinalize();
}

inline static FILE *wmFileOpen()
{
  const int len= strlen(GlOpts.filename)+ 8;
  char *fn= (char *) malloc(sizeof(char)* len);
  memset(fn, 0, sizeof(char)* len);
  memcpy(fn, GlOpts.filename, sizeof(char)* (len- 8));
  /* printf("fn0= \"%s\"\n", fn); */
  strcat(fn, ".wm.dat");
  /* printf("fn1= \"%s\", fn[-3]= '%c', fn[-2]= '%c', fn[-1]= '%c'\n", fn, fn[len- 3], fn[len- 2], fn[len- 1]); */

  FILE *fp= fopen(fn, "w");
  FREE_IF(fn);
  return fp;
}
void DoMonitoringWindowed()
{
  const int np= TraceGetNumProcs();

  struct {
    double *since;              /* elapsed */
    double *crit;
    int *state;
  } last= { NULL, NULL };

  struct {
    struct {
      double *critical;
      double *useful;
    } prev;

    struct {
      double *critical;
      double *useful;
    } curr;

    double *nevts;

    double tMin;
    double tMax;
    double step;
  } bin= { { NULL, NULL }, { NULL, NULL }, NULL, 0.0, 0.0, 0.0 };

  last.since= (double *) malloc(sizeof(double)* np);
  memset(last.since, 0, sizeof(double)* np);

  last.crit= (double *) malloc(sizeof(double)* np);
  memset(last.crit, 0, sizeof(double)* np);

  last.state= (int *) malloc(sizeof(int)* np);
  memset(last.state, 0, sizeof(int)* np);

  bin.prev.critical= (double *) malloc(sizeof(double)* np);
  memset(bin.prev.critical, 0, sizeof(double)* np);

  bin.prev.useful= (double *) malloc(sizeof(double)* np);
  memset(bin.prev.useful, 0, sizeof(double)* np);

  bin.curr.critical= (double *) malloc(sizeof(double)* np);
  memset(bin.curr.critical, 0, sizeof(double)* np);

  bin.curr.useful= (double *) malloc(sizeof(double)* np);
  memset(bin.curr.useful, 0, sizeof(double)* np);

  bin.nevts= (double *) malloc(sizeof(double)* np);
  memset(bin.nevts, 0, sizeof(double)* np);

  const double t0= TraceGetProgStartTimeMin();
  const double t1= TraceGetProgEndTimeMax();

  for(int ip= 0; ip< np; ++ip) {
    bin.curr.critical[ip]= bin.curr.useful[ip]= t0;
    last.since[ip]= last.crit[ip]= t0; last.state[ip]= -1;
  }

  FILE *fp= wmFileOpen();
  fprintf(fp, "%15s %15s %15s %15s %15s %15s %15s %15s\n",
          "#elapsed-1", "max-ideal-2", "avg-ideal-3", "max-useful-4", "avg-useful-5",
          "elapsed-loc-6", "ideal-loc-7", "min-nevts-8");

  bin.step= GlOpts.win_mon.win_len;
  bin.tMin= t0;
  TraceResetProcIters();
  const double pfactor= 1.0/ ((double) np);
  const double nevts_threshold= sqrt((double) np); /* MIN(128.0,((double) np)); */
  const long nbins= (long) ceil((t1- t0)/ bin.step);
  for(long ibin= 0; ibin< nbins; ++ibin) {
    bin.tMax= MIN(bin.tMin+bin.step, t1);
    memcpy(bin.prev.critical, bin.curr.critical, sizeof(double)* np);
    memcpy(bin.prev.useful, bin.curr.useful, sizeof(double)* np);
    memset(bin.nevts, 0, sizeof(double)* np);

again:
    /* work till elapsed exceeds end-of-the-bin */
    for(int ip= 0; ip< np; ++ip) {
      for(; TraceRemainsProcEvts(ip); TraceIncrIterProcEvts(ip)) {
        const double t= TraceGetAtCurrProcEvt(ip);
        if(t> bin.tMax) {
          break;
        }

        const double delta= t- last.since[ip];
        if(0== last.state[ip]) { /* leaving useful */
          bin.curr.critical[ip]+= delta;
          bin.curr.useful[ip]+= delta;
        } else {                /* leaving MPI/special */
          if(TraceGetCritCurrProcEvt(ip)> 0.1) {
            bin.curr.critical[ip]+= TraceGetCritCurrProcEvt(ip)- last.crit[ip];
          } else {
            Error("%d: at %.0lf, critical is %.0lf (%d -> %d)\n", ip,
                  t, TraceGetCritCurrProcEvt(ip),
                  last.state[ip], TraceGetIdCurrProcEvt(ip));
          }
        }
        bin.nevts[ip]+= 1.0;

        last.since[ip]= t;
        if(TraceGetCritCurrProcEvt(ip)> 0.1) {
          last.crit[ip]= TraceGetCritCurrProcEvt(ip);
        }
        last.state[ip]= TraceGetIdCurrProcEvt(ip);
      }
    }

    double nevtsmin= DBL_MAX;
    for(int ip= 0; ip< np; ++ip) {
      nevtsmin= MIN(bin.nevts[ip], nevtsmin);
    }
#if 1
    if(nevtsmin< nevts_threshold&& bin.tMax< t1) {
      ++ibin;
      bin.tMax= MIN(bin.tMax+bin.step, t1);
      goto again;
    }
#endif

    /* fill the remaing gap till end-of-the-bin */
    double tcmax= 0.0, tcmin= DBL_MAX;
    for(int ip= 0; ip< np; ++ip) {
      const double tremains= bin.tMax- last.since[ip];
      switch(last.state[ip]) {
      case 0:
        bin.curr.critical[ip]+= tremains;
        bin.curr.useful[ip]+= tremains;
        break;
      default: {
          if(TraceGetCritCurrProcEvt(ip)> 0.1) {
            const double delta= TraceGetCritCurrProcEvt(ip)- last.crit[ip];
            if(delta> tremains) {
              bin.curr.critical[ip]+= tremains;
              last.crit[ip]+= tremains;
            } else {
              bin.curr.critical[ip]+= delta;
              last.crit[ip]+= delta;
            }
          } else {
          }
        }
        break;
      }

      last.since[ip]= bin.tMax;
      tcmax= MAX(bin.curr.critical[ip], tcmax);
      tcmin= MIN(bin.curr.critical[ip], tcmin);
    }

    double uavg= 0.0, umax= 0.0, cavg= 0.0, cmax= 0.0;
    for(int ip= 0; ip< np; ++ip) {
      const double u= bin.curr.useful[ip]- bin.prev.useful[ip];
      uavg+= u;
      umax= MAX(u, umax);

      const double c= bin.curr.critical[ip]- bin.prev.critical[ip];
      cavg+= c;
      cmax= MAX(c, cmax);
    }
    uavg*= pfactor;
    cavg*= pfactor;

    double crit= cmax;
    if(crit< umax) {
      Debug1("Forced correction: critical < max-useful at bin [%.0lf, %.0lf)\n",
             bin.tMin, bin.tMax);
      crit= umax;
    }

    if(crit> bin.tMax- bin.tMin) {
      Debug1("Forced correction: critical > elapsed at bin [%.0lf, %.0lf)\n",
             bin.tMin, bin.tMax);
      crit= bin.tMax- bin.tMin;
    }

    fprintf(fp, "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
            bin.tMax, cmax, cavg, umax, uavg, bin.tMax- bin.tMin, crit, nevtsmin);

    bin.tMin= bin.tMax;
  }
  fclose(fp); fp= NULL;

  FREE_IF(bin.nevts);
  FREE_IF(bin.curr.useful);
  FREE_IF(bin.curr.critical);
  FREE_IF(bin.prev.useful);
  FREE_IF(bin.prev.critical);
  FREE_IF(last.state);
  FREE_IF(last.crit);
  FREE_IF(last.since);
}
