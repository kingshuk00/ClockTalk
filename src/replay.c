/*
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#include"build_info.h"
#include"common.h"
#include"utils.h"
#include"arg_opt_parser.h"
#include"paraver_file.h"
#include"trace_data.h"
#include"paraver.h"
#include"clocks.h"
#include"collectives.h"
#include"monitoring.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>
#include<float.h>

inline static void initialiseClocks(const int np)
{
  ClockInit(np);
  ClockStart(np, TraceGetProgStartTimeMin(), TraceGetPtrProcTimeline());
}
static void initialiseCollectives(const int np)
{
  CollsAlloc();
  CollsResetAll();
}
inline static int checkEvtsCompletion(const int p, bool *const completed)
{
  completed[p]= !TraceRemainsProcEvts(p);
  return completed[p]? 1: 0;
}

/* aliases to avoid calling functions with long names */
inline static const char *evt2name(const int evtId) { return GetParaverMPIEvtName(evtId); }
inline static int cevt(const int p) { return TraceGetIdCurrProcEvt(p); }
inline static int pevt(const int p) { return TraceGetIdPrevProcEvt(p); }
inline static const char *cevtname(const int p) { return evt2name(cevt(p)); }
inline static const char *pevtname(const int p) { return evt2name(pevt(p)); }
inline static double tcevt(const int p) { return TraceGetAtCurrProcEvt(p); }
inline static double tpevt(const int p) { return TraceGetAtPrevProcEvt(p); }

inline static bool excuse(const int prev, const int curr)
{
  if(GlOpts.sim_opts.ignore.trace_evts&& (-4== prev|| -4== curr)) {
    return true;
  }

  if(GlOpts.sim_opts.ignore.flush_evts&& (-3== prev|| -3== curr)) {
    return true;
  }

  if(GlOpts.sim_opts.ignore.disabled_tracing&& ((-2== prev&& 0== curr)||
                                                -2== curr)) {
    return true;
  }

  ErrorIf(-99==curr, "Invalid state encountered\n");

  return false;
}

static void enterMPI_Init(const int p)
{
  ClockPauseMPI(p, tcevt(p), cevt(p));
  Debug1("%d: enter %s at %.0lf -> %.0lf\n", p, cevtname(p), tcevt(p),
         ClockGetCritical(p));
  TraceSetCritCurrProcEvt(p, ClockGetCritical(p));
  TraceIncrIterProcEvts(p);
}
static void leaveMPI_InitAt(const int p, const double at)
{
  ClockSetCritical(p, at);
  ClockPlay(p, TraceGetAtCurrProcEvt(p), 0);
  Debug1("%d: exit %s at %.0lf -> %.0lf\n", p, pevtname(p), tcevt(p), at);
  TraceSetCritCurrProcEvt(p, ClockGetCritical(p));
  TraceIncrIterProcEvts(p);
}

static void localiseMPIUntilEvtExcl(const int p, const int untilEvt)
{
  const long untilIx= TraceSearchIterProcEvtId(p, untilEvt);
  while(TraceGetIterProcEvts(p)< untilIx) {
    const int eprev= pevt(p);
    const int e= cevt(p);
    const double t= tcevt(p);

    if(!excuse(eprev, e)) {
      if(e> 0) {
        ClockPauseMPI(p, t, e);
      } else if(e< 0) {
        ClockPauseTrace(p, t, e);
      } else {
        ClockPlay(p, t, e);
      }
      TraceSetCritCurrProcEvt(p, ClockGetCritical(p));
    }
    TraceIncrIterProcEvts(p);
  }
}
inline static void localiseMPIUntilEvtIncl(const int p, const int untilEvt)
{
  localiseMPIUntilEvtExcl(p, untilEvt);
  const long untilIx= TraceGetIterProcEvts(p)+ 2;
  while(TraceGetIterProcEvts(p)< untilIx) {
    TraceIncrIterProcEvts(p);
  }
}
static void playMPI_Init(const int np)
{
  if(TraceAllHaveMPIInitEvt()) {
    double lastEntry= 0.0;
    for(int ip= 0; ip< np; ++ip) {
      localiseMPIUntilEvtExcl(ip, 31);
      enterMPI_Init(ip);
#if 0
      ClockDebug(ip);
#endif
      lastEntry= MAX(ClockGetCritical(ip), lastEntry);
    }
    for(int ip= 0; ip< np; ++ip) {
      leaveMPI_InitAt(ip, lastEntry);
    }
  } else if(TraceAnyHasMPIInitEvt()) {
    for(int ip= 0; ip< np; ++ip) {
      if(TraceProcHasMPIInitEvt(ip)) {
        localiseMPIUntilEvtIncl(ip, 31);
      }
    }
  }
}

static int nstucks= 0;
inline static bool posted(const double *const t) { return t[2]> 0.1; }
inline static bool settled(const double *const t) { return t[2]< -0.1; }
inline static bool seen(const double *const t) { return posted(t)|| settled(t); }
inline static const char *seenstatus(const double *const t) { return (posted(t)? "posted": "settled"); }
inline static bool instant(const double *const t) { return SameTime(t[1], t[0]); }
inline static int sremote(const int p, const long ix) { return TraceGetProcSendRemote(p, ix); }
inline static int rremote(const int p, const long ix) { return TraceGetProcRecvRemote(p, ix); }
inline static bool rillogical(const int p, const long ix) { return TraceGetSendAtProcRecv(p, ix, 0)> TraceGetProcRecvAt(p, ix, 1); }
static void postOneSend(const int p, const long ix)
{
  double *const tsends= TraceGetProcSendAts(p, ix);
  if(seen(tsends)) {
    Debug1("%d: send->%d: event at %.0lf, already %s at %.0lf.\n", p, sremote(p,
           ix),
           tcevt(p), seenstatus(tsends), fabs(tsends[2]));
    ErrorIf(!SameTime(fabs(tsends[2]), ClockGetCritical(p)),
            "%d: send->%d: post override %.0lf->%.0lf\n", p, sremote(p, ix),
            fabs(tsends[2]), ClockGetCritical(p));
  }
  ErrorIf(tcevt(p)!= tsends[0], "%d: evt-time %.0lf but send-time %.0lf\n", p,
          tcevt(p), tsends[0]);
  Debug1("%d: send-%d: start %.0lf -> %.0lf\n", p, sremote(p, ix), tcevt(p),
         ClockGetCritical(p));
  TraceSetProcSendAt(p, ix, 2, ClockGetCritical(p));
}
static void postSends(const int p)
{
  IndexList *ixs= TraceGetCurrProcEvtSends(p, 0);
  while(NULL!= ixs) {
    postOneSend(p, ixs->i);
    ixs= ixs->next;
  }
}
static void postOneRecv(const int p, const long ix)
{
  double *const trecvs= TraceGetProcRecvAts(p, ix);
  if(seen(trecvs)) {
    Debug1("%d: recv<-%d: event at %.0lf, already %s at %.0lf.\n", p, rremote(p,
           ix),
           tcevt(p), seenstatus(trecvs), fabs(trecvs[2]));
    ErrorIf(!SameTime(fabs(trecvs[2]), ClockGetCritical(p)),
            "%d: recv<-%d: post override %.0lf->%.0lf\n", p, rremote(p, ix),
            fabs(trecvs[2]), ClockGetCritical(p));
  }
  ErrorIf(tcevt(p)!= trecvs[0], "%d: evt-time %.0lf but recv-time %.0lf\n", p,
          tcevt(p), trecvs[0]);
  Debug1("%d: recv-%d: start %.0lf -> %.0lf\n", p, rremote(p, ix), tcevt(p),
         ClockGetCritical(p));
  TraceSetProcRecvAt(p, ix, 2, ClockGetCritical(p));
}
static void postRecvs(const int p)
{
  IndexList *ixr= TraceGetCurrProcEvtRecvs(p, 0);
  while(NULL!= ixr) {
    postOneRecv(p, ixr->i);
    ixr= ixr->next;
  }
}
static void postMsgs(const int p)
{
  postSends(p);
  postRecvs(p);
}
static int postColls(const int p)
{
  int ret= 0;
  if(!ParaverCollEvtIsDimemasCompliant(cevt(p))) {
    goto bye;
  }

  const int c= TraceGetCommCurrProcColl(p);
  if(TraceIsCommSelf(c)) {
    goto bye;
  }

  const double t= tcevt(p);
  if(IsCollAvailable(c, p)) {
    Debug1("%d: coll %s(%d): start %.0lf -> unresolved resource-unavailable\n", p,
           cevtname(p), c, t);
    goto noresolve;
  }

  ErrorIf(fabs(TraceGetAtCurrProcColl(p, 0)- t)> 0.1,
          "%d: coll-evt at %.0lf: conflicting coll-record at %.0lf\n", p, t,
          TraceGetAtCurrProcColl(p, 0));

  EnterColl(c, p, t, cevt(p));
  goto bye;

noresolve:
  ret= 1;
bye:
  return ret;
}
static int settleColls(const int p)
{
  if(!ParaverCollEvtIsDimemasCompliant(pevt(p))) {
    return 0;
  }
  int ret= 0;
  if(!TraceIsCommSelf(TraceGetCommCurrProcColl(p))) {
    ret= LeaveColl(TraceGetCommCurrProcColl(p), p, tcevt(p), pevt(p));
  }

  if(0== ret) {
    TraceIncrIterProcColls(p);
  }
  return ret;
}

inline static bool NonblockingSendExit(const int p) { return 3== pevt(p)|| 36== pevt(p)|| 37== pevt(p)|| 38== pevt(p); }
inline static bool classicNonblockingSendExit(const int p,
                                              const double *const t) { return NonblockingSendExit(p)&& SameTime(t[0], tpevt(p))&& SameTime(t[1], tcevt(p)); }

/* returns 1 if this send is settled, 0 otherwise */
static int settleOneSend(const int p, const long ix)
{
  int ret= 1; const char *msg= "default";
  double *const tsends= TraceGetProcSendAts(p, ix);
  ErrorIf(tpevt(p)> tsends[1]&& tcevt(p)< tsends[1],
          "%d: evt-time %.0lf:%.0lf but send-time %.0lf\n", p, tpevt(p), tcevt(p),
          tsends[1]);

  if(settled(tsends)) {
    goto bye;  /* already settled: just return */
  }

  if(instant(tsends)) {         /* instant in trace: settle */
    msg= "instant"; goto settle;
  }

  if(classicNonblockingSendExit(p, tsends)) {     /* Isend-exit: settle */
    msg= "non-blocking"; goto settle;
  }

  if(TraceGetProcSendSize(p, ix)< GlOpts.sim_opts.eager_limit) {
    msg= "eager"; goto settle;
  }

  const double trecv= fabs(TraceGetRecvAtProcSend(p, ix, 2));
  if(trecv> 0.1) {           /* remote recv is posted, settle */
    if(trecv> ClockGetCritical(
         p)) { /* recv posted after clock-val - serialisation */
      ClockSetCritical(p, trecv);
    }
    msg= "rendezvous"; goto settle;
  }

  Debug1("%d: send-%d: end %.0lf -> unresolved rendezvous\n", p, sremote(p, ix),
         tcevt(p));
  ret= 0; goto bye;

settle:
  TraceSetProcSendAt(p, ix, 2, -tsends[2]); /* flip sign */
  Debug1("%d: send-%d: end %.0lf -> %.0lf %s\n", p, sremote(p, ix), tcevt(p),
         ClockGetCritical(p), msg);

bye:
  return ret;
}
/* returns 0 if all sends concluding at this event are settled */
/* returns 1 otherwise: due to corresponding recvs not yet posted */
static int settleSends(const int p)
{
  IndexList *ixs= TraceGetCurrProcEvtSends(p, 1);
  int total= 0, settled= 0;
  while(NULL!= ixs) {
    ++total;
    settled+= settleOneSend(p, ixs->i);
    ixs= ixs->next;
  }
  return (total== settled? 0: 1);
}

static int settleOneRecv(const int p, const long ix)
{
  int ret= 1; const char *msg= "default";
  double *const trecvs= TraceGetProcRecvAts(p, ix);
  ErrorIf(tpevt(p)> trecvs[1]&& tcevt(p)< trecvs[1],
          "%d: evt-time %.0lf:%.0lf but recv-time %.0lf\n", p, tpevt(p), tcevt(p),
          trecvs[1]);

  if(settled(trecvs)) {
    goto bye;  /* already settled */
  }

  /* if(instant(trecvs)) { msg= "instant";  goto settle; } */

  const double tsend= fabs(TraceGetSendAtProcRecv(p, ix, 2));
  if(tsend> 0.1) {           /* remote send is posted, settle */
    if(tsend> ClockGetCritical(
         p)) { /* send posted after clock-val - serialisation */
      ClockSetCritical(p, tsend);
    }
    msg= "remote-post"; goto settle;
  }

  /* If all are stuck and unresolved illogical receives are encountered,
   * settle them with current time.
   * Illogical receives finish before corresponding sends are even posted.
   * They arise due to clock-skew and/or time dilation and are usually fixable.
   * Rarely, incorrect message matching in Paraver files show up this way,
   * at times crossing other causally linked events. They are not fixable.
   * This "hack" ends such receives so that the best possible monitor even with
   * a mismatched trace can be obtained.
   */
  if(nstucks> 0&& rillogical(p, ix)) {
    Error("%d: illogical send-recv pair (%.0lf:%.0lf) <- rank-%d (%.0lf:%.0lf)\n",
          p,
          trecvs[0], trecvs[1], rremote(p, ix),
          TraceGetSendAtProcRecv(p, ix, 0), TraceGetSendAtProcRecv(p, ix, 1));
    msg= "halted illogical recv"; goto settle;
  }

  Debug1("%d: recv-%d: end %.0lf -> unresolved remote-not-posted\n", p, rremote(p,
         ix), tcevt(p));
  ret= 0; goto bye;

settle:
  TraceSetProcRecvAt(p, ix, 2, -ClockGetCritical(p));
  Debug1("%d: recv-%d: end %.0lf -> %.0lf %s\n", p, rremote(p, ix), tcevt(p),
         ClockGetCritical(p), msg);

bye:
  return ret;
}
/* returns 0 if all recvs concluding at this event are settled */
/* returns 1 otherwise: due to corresponding sends not yet posted */
static int settleRecvs(const int p)
{
  IndexList *ixr= TraceGetCurrProcEvtRecvs(p, 1);
  int total= 0, settled= 0;
  while(NULL!= ixr) {
    ++total;
    settled+= settleOneRecv(p, ixr->i);
    ixr= ixr->next;
  }
  return (total== settled? 0: 1);
}
static int settleMsgs(const int p)
{
  if(0!= settleSends(p)) {
    return 1;
  }
  if(0!= settleRecvs(p)) {
    return 1;
  }
  return 0;
}

/* progresses as much as possible without talking */
#if 1
static int processRank(const int p)
{
  int movement= 0;
  while(TraceRemainsProcEvts(p)) {
    const double t= tcevt(p);
    const int e= cevt(p);

    if(excuse(pevt(p), e)) {
      goto excused;
    }

    if(e> 0) {                       /* enter MPI */
      ClockPauseMPI(p, t, e);
      postMsgs(p);
      if(0!= settleMsgs(p)) {
        Debug1("%d: p2p: all msgs are not settled\n", p);
        break;
      }
      if(0!= postColls(p)) {
        Debug1("%d: coll: not posted - another collective on comm\n", p);
        break;
      }
    } else if(e< 0) {                /* enter trace-events */
      ClockPauseTrace(p, t, e);
    } else {                         /* enter useful */
      if(0!= settleColls(p)) {
        Debug1("%d: coll: not settled - wait for others\n", p);
        break;
      }
      postMsgs(p);
      if(0!= settleMsgs(p)) {
        Debug1("%d: p2p: all msgs are not settled\n", p);
        break;
      }
      ClockPlay(p, t, e);
    }
    TraceSetCritCurrProcEvt(p, ClockGetCritical(p));

    const char *status= "processed";
    goto through;

excused:
    /* TraceSetCurrEvtCrit(z); */
    status= "excused";

through:
    Debug1("%d: %s -> %s at %.0lf (%.0lf) - %s\n", p, pevtname(p),
           cevtname(p), ClockGetElapsed(p), ClockGetCritical(p), status);
    TraceIncrIterProcEvts(p);
    ++movement;
  }
  return movement;
}
#endif

/* progresses one step */
#if 0
static int processRank(const int p)
{
  if(!TraceRemainsProcEvts(p)) {
    return 0;
  }
  int movement= 0;

  const double t= tcevt(p);
  const int e= cevt(p);

  if(unlikely(excuse(pevt(p), e))) {
    goto excused;
  }

  if(e> 0) {                       /* enter MPI */
    ClockPauseMPI(p, t, e);
    postMsgs(p);
    if(0!= settleMsgs(p)) {
      Debug1("%d: p2p: all msgs are not settled\n", p);
      return movement;
    }
    if(0!= postColls(p)) {
      Debug1("%d: coll: not posted - another collective on comm\n", p);
      return movement;
    }
  } else if(e< 0) {                /* enter trace-events */
    ClockPauseTrace(p, t, e);
  } else {                         /* enter useful */
    if(0!= settleColls(p)) {
      Debug1("%d: coll: not settled - wait for others\n", p);
      return movement;
    }
    postMsgs(p);
    if(0!= settleMsgs(p)) {
      Debug1("%d: p2p: all msgs are not settled\n", p);
      return movement;
    }
    ClockPlay(p, t, e);
  }

  const char *status= "processed";
  goto through;

excused:
  status= "excused";

through:
  Debug1("%d: %s -> %s at %.0lf (%.0lf) - %s\n", p, pevtname(p),
         cevtname(p), ClockGetElapsed(p), ClockGetCritical(p), status);
  TraceIncrIterProcEvts(p);
  ++movement;

  return movement;
}
#endif

#if 0
/* monitoring rank's timepoints */
static void calcMonRanksTimepoints(long *const n, double *const umax,
                                   double *uavg)
{
  const int p= GlOpts.evt_mon.rank;
  long npoints= 0;

  TraceResetProcIters();
  while(TraceRemainsProcEvts(p)) {
    const int e= cevt(p);
    if(excuse(pevt(p), e)) {
      goto excused;
    }
    if(0== e) {
      ++npoints;
    }

    goto through;

excused:
    ;

through:
    TraceIncrIterProcEvts(p);
  }
  printf("#points in tl-file: %ld\n", npoints);
  TraceResetProcIters();
}
#endif

static void processTrace()
{
  /* calcMonRanksTimepoints(NULL, NULL, NULL); */
  const double t0= Timer_s();
  TraceConnectEvtsToMsgs();
  if(GlOpts.show_opts.timings) {
    printf("Connecting MPI events to p2p calls took %.1lf s\n", Timer_s()- t0);
  }

  const int np= TraceGetNumProcs();

  initialiseClocks(np);
  initialiseCollectives(np);

  TraceResetProcIters();

  playMPI_Init(np);

  int ncompleted= 0;
  bool *completed= (bool *) malloc(sizeof(bool)* np);
  for(int ip= 0; ip< np; ++ip) {
    ncompleted+= checkEvtsCompletion(ip, completed);
  }

  long *ix= (long *) malloc(sizeof(long)* np);
  memset(ix, 0, sizeof(long)* np);

  const int maxstucks= np;
  while(ncompleted< np) {       /* main loop over events records */
    int movement= 0;
    for(int ip= 0; ip< np; ++ip) {
      if(completed[ip]) {
        Debug1("%d: events completed, skipping\n", ip);
        continue;
      }

      movement+= processRank(ip);

      ncompleted+= checkEvtsCompletion(ip, completed);
    }

    if(0== movement&& false) {
      Debug1("This iteration has not progressed\n");
      ++nstucks;
      if(nstucks> maxstucks- 1) {
        Error("All are stuck for %d iterations, getting out\n", maxstucks);
        break;
      }
    } else {
      nstucks= 0;
    }
  }

  if(maxstucks!= nstucks) {
#if 0
    bool appEnds= false;        /* 40000001:0 doesn't exist */
    for(int ip= 0; ip< np; ++ip) {
      if(tpevt(ip)< TraceGetProcEndTime(ip)) { /* 40000001:0 exists */
        appEnds= true;
        break;
      }
    }
#endif

    const double tEnd= TraceGetProgEndTimeMax();
    /* if(appEnds) { */
    if(true) {
      for(int ip= 0; ip< np; ++ip) {
        ClockEnd(ip, TraceGetProcEndTime(ip), tEnd);
      }
    } else {
      for(int ip= 0; ip< np; ++ip) {
        ClockEnd(ip, tEnd, tEnd);
      }
    }
  } else {
    Error("Wrong results!!\n");
  }
}

static void showStats()
{
  const int np= TraceGetNumProcs();
  const double n2u= 1.0e-3;
  const double runtime= ClockGetMaxElapsed(np)* n2u;
  const double runtime_inv= 1.0/ runtime;
  const double runtime_traced= ClockGetMaxTraced(np)* n2u;
  const double runtime_traced_inv= 1.0/ runtime_traced;
  const double runtime_traced_ideal= (ClockGetMaxCritical(np))* n2u;
  const double useful_max= ClockGetMaxUseful(np)* n2u;
  const double useful_avg= ClockGetAvgUseful(np)* n2u;

  FILE *fp= stdout;
  if(GlOpts.show_opts.pretty) {
    fprintf(fp, "==============================================\n");
    fprintf(fp, "       runtime= %.2lf us\n", runtime);
    fprintf(fp, "traced-runtime= %.2lf us\n", runtime_traced);
    fprintf(fp, " ideal-runtime= %.2lf us\n", runtime_traced_ideal);

#if 0
    const double diff= runtime_traced_ideal- 421272658.40; /* from rawdata.csv */
    if(fabs(diff)> 0.1) {
      fprintf(fp, " (%.2lf, %.2lf%%)", diff, fabs(diff)* 100.0/ runtime_traced_ideal);
    } else {
      fprintf(fp, " (match)");
    }
#endif

    fprintf(fp, "    max-useful= %.2lf us\n", useful_max);
    fprintf(fp, "    avg-useful= %.2lf us\n", useful_avg);
    fprintf(fp, "==============================================\n");
    fprintf(fp, "Overview of the Efficiency metrics:\n");
    fprintf(fp, "==============================================\n");
    fprintf(fp, "                       Trace mode |        MPI\n");
    fprintf(fp, "          Processes [Trace Order] |  %6d[1]\n", np);
    fprintf(fp, "==============================================\n");
    fprintf(fp, "Global efficiency                 |    %6.2lf%% (%3.0lf%%)\n",
            useful_avg* 100.0* runtime_inv, useful_avg* 100.0* runtime_traced_inv);
    fprintf(fp, "-- Parallel efficiency            |    %6.2lf%% (%3.0lf%%)\n",
            useful_avg* 100.0* runtime_inv, useful_avg* 100.0* runtime_traced_inv);
    fprintf(fp, "   -- Load balance                |    %6.2lf%%\n",
            useful_avg* 100.0/ useful_max);
    fprintf(fp, "   -- Communication efficiency    |    %6.2lf%% (%3.0lf%%)\n",
            useful_max* 100.0* runtime_inv, useful_max* 100.0* runtime_traced_inv);
    fprintf(fp, "      -- Serialization efficiency |    %6.2lf%%\n",
            useful_max* 100.0/ runtime_traced_ideal);
    fprintf(fp, "      -- Transfer efficiency      |    %6.2lf%% (%3.0lf%%)\n",
            runtime_traced_ideal* 100.0* runtime_inv,
            runtime_traced_ideal* 100.0* runtime_traced_inv);
    fprintf(fp, "-- Computation scalability        |  Non-Avail\n");
    fprintf(fp, "   -- IPC scalability             |  Non-Avail\n");
    fprintf(fp, "   -- Instruction scalability     |  Non-Avail\n");
    fprintf(fp, "   -- Frequency scalability       |  Non-Avail\n");
    fprintf(fp, "==============================================\n");
  } else {
    fprintf(fp,
            "runtime= %.2lf us\n traced= %.2lf us\n  ideal= %.2lf us\n   uavg= %.2lf us\n   umax= %.2lf us\n",
            runtime, runtime_traced, runtime_traced_ideal, useful_avg, useful_max);
    fprintf(fp, "pe= %.4lf, lb= %.4lf, ser= %.4lf, trf= %.4lf\n",
            useful_avg* runtime_inv, useful_avg/ useful_max,
            useful_max/ runtime_traced_ideal, runtime_traced_ideal* runtime_inv);
  }
}

inline static void PrintGlobalOpts()
{
#if 0
  printf("GlOpts:\n  filename: \"%s\"\n\n", GlOpts.filename);
  printf("  show_opts:\n");
  printf("    diag: %d\n", GlOpts.show_opts.diag);
  printf("    error: %d\n", GlOpts.show_opts.error);
  printf("    io_timings: %s\n", GlOpts.show_opts.timings? "true": "false");
  printf("    profile: %s\n", GlOpts.show_opts.profile? "true": "false");
  printf("    pretty: %s\n", GlOpts.show_opts.pretty? "true": "false");
  printf("\n  win_mon:\n");
  printf("    win_len: %.6e\n", GlOpts.win_mon.win_len);
  printf("    nwins_sma: %d\n", GlOpts.win_mon.nwins_sma);
  printf("    enabled: %s\n", GlOpts.win_mon.enabled? "true": "false");
  printf("  \n  evt_mon:\n");
  printf("    rank: %d\n", GlOpts.evt_mon.rank);
  printf("    nevts_report: %d\n", GlOpts.evt_mon.nevts_report);
  printf("    enabled: %s\n", GlOpts.evt_mon.enabled? "true": "false");
  printf("  \n  sim_opts:\n");
  printf("    eager_limit: %.0lf\n", GlOpts.sim_opts.eager_limit);
  printf("    ignore:\n");
  printf("      trace_evts: %s\n",
         GlOpts.sim_opts.ignore.trace_evts? "true": "false");
  printf("      flush_evts: %s\n",
         GlOpts.sim_opts.ignore.flush_evts? "true": "false");
  printf("      disabled_tracing: %s\n",
         GlOpts.sim_opts.ignore.disabled_tracing? "true": "false");

  /* exit(0); */
#endif
}

#if 0
void PrintAbc(const int ip)
{
  printf("%.0lf -> %.0lf -> %.0lf\n", tpevt(ip), tcevt(ip), tnevt(ip));
  printf("%.0lf -> %.0lf -> %.0lf\n\n", ctpevt(ip), ctcevt(ip), ctnevt(ip));
  fflush(stdout);
}
#endif

#if 0
void Abc()
{
  const int np= TraceGetNumProcs();
  struct {
    double *tprv;
    double *tsim;
    int *evt;
    double *at;
    bool *isEnabled;
  } last= { NULL, NULL, NULL, NULL, NULL };

  struct {
    double tMin;
    double tMax;
    double step;
    double cMaxLast;
    double *nevts;
    double *useful;
    double *critic;
  } bin= { 0.0, 0.0, 0.0, 0.0, NULL, NULL, NULL };

  last.tprv= (double *) malloc(sizeof(double)* np);
  memset(last.tprv, 0, sizeof(double)* np);

  last.tsim= (double *) malloc(sizeof(double)* np);
  memset(last.tsim, 0, sizeof(double)* np);

  last.evt= (int *) malloc(sizeof(int)* np);
  memset(last.evt, 0, sizeof(int)* np);

  last.at= (double *) malloc(sizeof(double)* np);
  memset(last.at, 0, sizeof(double)* np);

  last.isEnabled= (bool *) malloc(sizeof(bool)* np);
  memset(last.isEnabled, 0, sizeof(bool)* np);

  bin.nevts= (double *) malloc(sizeof(double)* np);
  memset(bin.nevts, 0, sizeof(double)* np);

  bin.useful= (double *) malloc(sizeof(double)* np);
  memset(bin.useful, 0, sizeof(double)* np);

  bin.critic= (double *) malloc(sizeof(double)* np);
  memset(bin.critic, 0, sizeof(double)* np);

  const double t0= TraceGetProgStartTimeMin();
  const double t1= TraceGetProgEndTimeMax();

  for(int ip= 0; ip< np; ++ip) {
    last.tprv[ip]= last.tsim[ip]= t0;
    last.evt[ip]= -1;
  }

  FILE *fp= fopen("abc.dat", "w");
  fprintf(fp, "#%14s %15s %15s %15s %15s %15s %15s %15s\n",
          "t1-1", "uavg-2", "umax-3", "cavg-4", "cmax-5", "elapsed-6", "crit-7",
          "min-nevts-8");

  bin.step= 2.0e9;
  bin.tMin= t0;
  bin.cMaxLast= t0;
  TraceResetProcIters();
  const double pfactor= 1.0/ ((double) np);
  const double nevts_threshold= sqrt((double) np); /* MIN(128.0,((double) np)); */
  const long nbins= (long) ceil((t1- t0)/ bin.step);
  for(long ibin= 0; ibin< nbins; ++ibin) {
    memset(bin.nevts, 0, sizeof(double)* np);
    memset(bin.useful, 0, sizeof(double)* np);
    memset(bin.critic, 0, sizeof(double)* np);
    bin.tMax= MIN(bin.tMin+bin.step, t1);

again:
    /* work till elapsed exceeds end-of-the-bin */
    for(int ip= 0; ip< np; ++ip) {
      for(; TraceRemainsProcEvts(ip); TraceIncrIterProcEvts(ip)) {
        if(tcevt(ip)> bin.tMax) {
          break;
        }

        if(0== last.evt[ip]) {
          bin.useful[ip]+= tcevt(ip)- last.tprv[ip];
        }
        bin.critic[ip]+= ctcevt(ip)- last.tsim[ip];
        bin.nevts[ip]+= 1.0;

        last.tprv[ip]= tcevt(ip);
        last.tsim[ip]= ctcevt(ip);
        last.evt[ip]= cevt(ip);
      }
    }

    double nevtsmin= DBL_MAX;
    for(int ip= 0; ip< np; ++ip) {
      nevtsmin= MIN(bin.nevts[ip], nevtsmin);
    }
#if 1
    if(nevtsmin< nevts_threshold&& bin.tMax< t1) {
      ++ibin;
      bin.tMax+= bin.step;
      goto again;
    }
#endif

    /* fill the remaing gap till end-of-the-bin */
    double tcmax= 0.0;
    for(int ip= 0; ip< np; ++ip) {
      const double tremains= bin.tMax- last.tprv[ip];
      switch(last.evt[ip]) {
      case 0:
        bin.useful[ip]+= tremains;
        bin.critic[ip]+= tremains;
        last.tsim[ip]+= tremains;
        break;
      default: {
          const double tcritnext= ctcevt(ip)- last.tsim[ip];
          if(tcritnext> tremains) {
            bin.critic[ip]+= tremains;
            last.tsim[ip]+= tremains;
          } else {
            bin.critic[ip]+= tcritnext;
            last.tsim[ip]+= tcritnext;
          }
        }
        break;
      }

      last.tprv[ip]= bin.tMax;
      tcmax= MAX(last.tsim[ip], tcmax);
    }

    double uavg= 0.0, umax= 0.0, cavg= 0.0, cmax= 0.0;
    for(int ip= 0; ip< np; ++ip) {
      uavg+= bin.useful[ip];
      umax= MAX(bin.useful[ip], umax);

      cavg+= bin.critic[ip];
      cmax= MAX(bin.critic[ip], cmax);
    }
    uavg*= pfactor;
    cavg*= pfactor;

    double crit= tcmax- bin.cMaxLast; bin.cMaxLast= tcmax;

    if(crit< umax) {
      printf("Forcing: critical was less than max-useful\n");
      crit= umax;
    }

    if(crit> bin.tMax- bin.tMin) {
      printf("Forcing: critical was more than time-window\n");
      crit= bin.tMax- bin.tMin;
    }

    fprintf(fp, "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
            bin.tMax, uavg, umax, cavg, cmax, bin.tMax- bin.tMin, crit, nevtsmin);

    bin.tMin= bin.tMax;
  }

  fclose(fp); fp= NULL;
}
#endif

int main(int argc, char *argv[])
{
  if(0!= ParseArgs(argc, argv)) {
    ArgHelp();
    return 0;
  }
  PrintGlobalOpts();

  Debug1("Running program built on %s at %s\n", CLOCKTALK_BUILD_DATE,
         CLOCKTALK_BUILD_TIME);

  const double t0= Timer_s();
  if(0!= ReadParaverFile(GlOpts.filename)) {
    Error("Problem reading paraver file \"%s\"\n", argv[1]);
    return 0;
  }
  const double t1= Timer_s();

  if(GlOpts.show_opts.timings) {
    printf("Reading Paraver file took %.1lf s\n", t1- t0);
  }

  processTrace();
  if(false) {
    FILE *fp= fopen("checking.txt", "w");
    for(TraceResetIterEvts(); TraceGetIterEvts()< TraceGetNumEvts();
        TraceIncrIterEvts()) {
      fprintf(fp, "%.9e %.9e %3d\n",
              TraceGetCurrEvtAt(), TraceGetCurrEvtCrit(), TraceGetCurrEvtId());
    }
    fclose(fp); fp= NULL;
  }

  showStats();

  ClockFinalize();

  if(GlOpts.show_opts.timings) {
    const double t2= Timer_s();
    printf("Replay took %.1lf s (total %.1lf s)\n", t2- t1, t2- t0);
  }

  if(GlOpts.evt_mon.enabled) {
    DoMonitoringEventBased();
  }

  if(GlOpts.win_mon.enabled) {
    DoMonitoringWindowed();
  }

  FREE_IF(GlOpts.filename);

  return 0;
}
