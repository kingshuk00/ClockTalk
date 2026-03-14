/*
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef REPLAY_TRACE_DATA_H__
#define REPLAY_TRACE_DATA_H__

#include<string.h>              /* memset() */
#include<math.h>                /* fabs */
#include<stdio.h>               /* printf() */
#include<stdlib.h>
#include<stdbool.h>

typedef struct IndexList_type__ {
  long i;
  struct IndexList_type__ *next;
} IndexList;
typedef struct ProcMap_struct__ {
  long *nums;                   /* len= #procs */
  long *iters;                  /* len= #procs */
  long **gids;
} ProcMap;
typedef struct TraceData_struct__ {
  long long runtime;
  char timeunit[4];
  int numnodes;                 /* #resources */
  int numapps;
  int numprocs;

  double extremities[2];

  struct {
    double (*extents)[2]; /* len= #procs; aka app-event 40000001:[1/0]  */
    double *tcomp;        /* len= #procs; aka 'useful' */
    double *tmpi;         /* len= #procs; out of MPI= 'runtime'- 'tmpi' */
    double *tflush;       /* len= #procs */
    double *disabledAt;   /* len= #procs */
    double *tdisabled;    /* len= #procs */
    bool *hasTraceInit;   /* len= #procs */
    bool *hasMPIInit;     /* len= #procs */
  } timeline;

  struct {
    long num;                   /* #communicators in whole trace */
    int *sizes;
    int **ranks;
  } comms;

  struct {
    long num;                   /* #events in whole trace */
    long iter;
    double *at;
    int *id; /* -1: end, -2: disabled, -3: flush, -4: trace-init, -99: invalid */
    double *crit;
    int *proc;

    IndexList **slist[2];
    IndexList **rlist[2];
  } evts;

  struct {
    long num;                   /* #communications in whole trace */
    long iter;
    double (*st)[3];
    double (*rt)[3];
    int *srank;
    int *rrank;
    double *size;
    int *tag;
    IndexList *ipsmem;
    IndexList *iprmem;
  } msgs;

  ProcMap pevts;
  ProcMap psends;
  ProcMap precvs;

  struct {
    long allnum;
    long *nums;
    long *iters;
    double (**at)[2];
    int **comm;
  } pcolls;

} TraceData;

struct ParaverFile_struct__;
extern TraceData *CreateTrace(const struct ParaverFile_struct__ *const);

extern TraceData *Trace0;
inline static void SetWorkingTrace(TraceData *const t) { Trace0= t; }

inline static long long TraceGetRuntime() { return Trace0->runtime; }
inline static const char *TraceGetTimeUnit() { return Trace0->timeunit; }
inline static int TraceGetNumNodes() { return Trace0->numnodes; }
inline static int TraceGetNumApps() { return Trace0->numapps; }
inline static int TraceGetNumProcs() { return Trace0->numprocs; }

/* timeline-proc-extents */
inline static double *TraceGetPtrProcTimeline() { return Trace0->timeline.extents[0]; }
inline static void TraceSetProcTimeline(const int p, const int i,
                                        const double t) { Trace0->timeline.extents[p][i]= t; }
inline static void TraceStartProcTimeline(const int p, const double at) { TraceSetProcTimeline(p, 0, at); }
inline static void TraceEndProcTimeline(const int p, const double at) { TraceSetProcTimeline(p, 1, at); }
inline static double TraceGetProcTimeline(const int p, const int i) { return Trace0->timeline.extents[p][i]; }
inline static double TraceGetProcStartTime(const int p) { return TraceGetProcTimeline(p, 0); }
inline static double TraceGetProcEndTime(const int p) { return TraceGetProcTimeline(p, 1); }

/* extremities */
inline static double *TraceGetPtrExtremities() { return Trace0->extremities; }
inline static void TraceCalculateTimelineExtremes()
{
  double *const ts= TraceGetPtrExtremities();
  ts[0]= TraceGetProcStartTime(0);
  ts[1]= TraceGetProcEndTime(0);
  for(int ip= 1; ip< TraceGetNumProcs(); ++ip) {
    if(TraceGetProcStartTime(ip)< ts[0]) {
      ts[0]= TraceGetProcStartTime(ip);
    }
    if(TraceGetProcEndTime(ip)> ts[1]) {
      ts[1]= TraceGetProcEndTime(ip);
    }
  }
}
inline static double TraceGetProgStartTimeMin() { return Trace0->extremities[0]; }
inline static double TraceGetProgEndTimeMax() { return Trace0->extremities[1]; }

/* proc-timeline-comp-duration */
inline static void TraceStartProcComp(const int p, const double at) { Trace0->timeline.tcomp[p]-= at; }
inline static void TraceEndProcComp(const int p, const double at) { Trace0->timeline.tcomp[p]+= at; }
inline static void TraceAddProcComp(const int p, const double t) { Trace0->timeline.tcomp[p]+= t; }
inline static double TraceGetProcCompDuration(const int p) { return Trace0->timeline.tcomp[p]; }

/* proc-timeline-mpi-duration */
inline static void TraceStartProcMPI(const int p, const double at) { Trace0->timeline.tmpi[p]-= at; }
inline static void TraceEndProcMPI(const int p, const double at) { Trace0->timeline.tmpi[p]+= at; }
inline static void TraceAddProcMPI(const int p, const double t) { Trace0->timeline.tmpi[p]+= t; }
inline static double TraceGetProcMPIDuration(const int p) { return Trace0->timeline.tmpi[p]; }

/* proc-timeline-flush-duraiton */
inline static void TraceStartProcFlush(const int p, const double at) { Trace0->timeline.tflush[p]-= at; }
inline static void TraceEndProcFlush(const int p, const double at) { Trace0->timeline.tflush[p]+= at; }
inline static double TraceGetProcFlushDuration(const int p) { return Trace0->timeline.tflush[p]; }

/* proc-timeline-disabled-duraiton */
inline static void TraceStartProcDisabled(const int p, const double at) { Trace0->timeline.disabledAt[p]= at; }
inline static void TraceEndProcDisabled(const int p, const double at)
{
  if(Trace0->timeline.disabledAt[p]> -0.1) {
    Trace0->timeline.tdisabled[p]+= at- Trace0->timeline.disabledAt[p];
    Trace0->timeline.disabledAt[p]= -1.0;
  }
}
inline static double TraceProcDisabledAt(const int p) { return Trace0->timeline.disabledAt[p]; } /* enabled: < 0; disabled: > 0 */
inline static double TraceGetProcDisabledDuration(const int p) { return Trace0->timeline.tdisabled[p]; }

/* proc-timeline-trace-init-is-present */
inline static void TraceSetTraceInitEvt(const int p) { Trace0->timeline.hasTraceInit[p]= true; }
inline static bool TraceHasTraceInitEvt(const int p) { return Trace0->timeline.hasTraceInit[p]; }

/* proc-timeline-MPI_Init-is-present */
inline static void TraceSetMPIInitEvt(const int p) { Trace0->timeline.hasTraceInit[p]= true; }
inline static bool TraceProcHasMPIInitEvt(const int p) { return Trace0->timeline.hasTraceInit[p]; }
inline static bool TraceAllHaveMPIInitEvt()
{
  for(int ip= 0; ip< Trace0->numprocs; ++ip) {
    if(!TraceProcHasMPIInitEvt(ip)) {
      return false;
    }
  }
  return true;
}
inline static bool TraceAnyHasMPIInitEvt()
{
  for(int ip= 0; ip< Trace0->numprocs; ++ip) {
    if(TraceProcHasMPIInitEvt(ip)) {
      return true;
    }
  }
  return false;
}

/* communicators */
inline static long TraceGetNumComms() { return Trace0->comms.num; }
inline static int *TraceGetPtrCommsSizes() { return Trace0->comms.sizes; }
inline static int *TraceGetPtrCommsRanks(const int ic) { return Trace0->comms.ranks[ic]; }
inline static int TraceGetCommRank(const int ic, const int i) { return Trace0->comms.ranks[ic][i]; }
inline static int TraceGetCommSize(const int ic) { return Trace0->comms.sizes[ic]; }
inline static bool TraceIsCommSelf(const int ic) { return 1== Trace0->comms.sizes[ic]; }
inline static int TraceGetSelfCommRank(const int ic) { return TraceGetCommRank(ic, 0); }

/* evts-num */
inline static void TraceSetNumEvts(const long num) { Trace0->evts.num= num; }
inline static long TraceGetNumEvts() { return Trace0->evts.num; }
/* evts-iter */
inline static void TraceResetIterEvts() { Trace0->evts.iter= 0; }
inline static void TraceIncrIterEvts() { ++(Trace0->evts.iter); }
inline static long TraceGetIterEvts() { return Trace0->evts.iter; }
inline static bool TraceRemainEvts() { return TraceGetIterEvts()< TraceGetNumEvts(); }
/* evts-at */
inline static void TraceSetPtrEvtsAt(double *at) { Trace0->evts.at= at; }
inline static double *TraceGetPtrEvtsAt() { return Trace0->evts.at; }
inline static void TraceSetEvtAt(const long it, const double t) { Trace0->evts.at[it]= t; }
inline static void TraceSetCurrEvtAt(const double t) { TraceSetEvtAt(TraceGetIterEvts(), t); }
inline static double TraceGetEvtAt(const long it) { return Trace0->evts.at[it]; }
inline static double TraceGetCurrEvtAt() { return TraceGetEvtAt(TraceGetIterEvts()); }
/* evts-id */
inline static void TraceSetPtrEvtsId(int *id) { Trace0->evts.id= id; }
inline static int *TraceGetPtrEvtsId() { return Trace0->evts.id; }
inline static void TraceSetEvtId(const long it, const int id) { Trace0->evts.id[it]= id; }
inline static void TraceSetCurrEvtId(const int id) { TraceSetEvtId(TraceGetIterEvts(), id); }
inline static int TraceGetEvtId(const long it) { return Trace0->evts.id[it]; }
inline static int TraceGetCurrEvtId() { return TraceGetEvtId(TraceGetIterEvts()); }
/* evts-crit */
inline static void TraceSetPtrEvtsCrit(double *crit) { Trace0->evts.crit= crit; }
inline static double *TraceGetPtrEvtsCrit() { return Trace0->evts.crit; }
inline static void TraceSetEvtCrit(const long it, const double crit) { Trace0->evts.crit[it]= crit; }
inline static void TraceSetCurrEvtCrit(const double crit) { TraceSetEvtCrit(TraceGetIterEvts(), crit); }
inline static double TraceGetEvtCrit(const long it) { return Trace0->evts.crit[it]; }
inline static double TraceGetCurrEvtCrit() { return TraceGetEvtCrit(TraceGetIterEvts()); }
/* evts-proc */
inline static void TraceSetPtrEvtsProc(int *proc) { Trace0->evts.proc= proc; }
inline static int *TraceGetPtrEvtsProc() { return Trace0->evts.proc; }
inline static void TraceSetEvtProc(const long it, const int proc) { Trace0->evts.proc[it]= proc; }
inline static void TraceSetCurrEvtProc(const int proc) { TraceSetEvtProc(TraceGetIterEvts(), proc); }
inline static int TraceGetEvtProc(const long it) { return Trace0->evts.proc[it]; }
inline static int TraceGetCurrEvtProc() { return TraceGetEvtProc(TraceGetIterEvts()); }

/* evts-slist */
inline static void TraceSetEvtSends(const int i, const long it,
                                    IndexList *sends) { Trace0->evts.slist[i][it]= sends; }
inline static IndexList **TraceGetPtrEvtSends(const int i) { return Trace0->evts.slist[i]; }
inline static IndexList *TraceGetEvtSends(const int i, const long it) { return Trace0->evts.slist[i][it]; }

/* evts-rlist */
inline static void TraceSetEvtRecvs(const int i, const long it,
                                    IndexList *recvs) { Trace0->evts.rlist[i][it]= recvs; }
inline static IndexList **TraceGetPtrEvtRecvs(const int i) { return Trace0->evts.rlist[i]; }
inline static IndexList *TraceGetEvtRecvs(const int i, const long it) { return Trace0->evts.rlist[i][it]; }

/* pevts-nums */
inline static long *TraceGetPtrNumProcEvts() { return Trace0->pevts.nums; }
inline static void TraceIncrNumProcEvts(const int p) { ++(Trace0->pevts.nums[p]); }
inline static long TraceGetNumProcEvts(const int p) { return Trace0->pevts.nums[p]; }
/* pevts-iters */
inline static long *TraceGetPtrIterProcEvts() { return Trace0->pevts.iters; }
inline static void TraceResetItersProcEvts() { memset(TraceGetPtrIterProcEvts(), 0, sizeof(long)* TraceGetNumProcs()); }
inline static void TraceSetIterProcEvts(const int p, const long it) { Trace0->pevts.iters[p]= it; }
inline static void TraceIncrIterProcEvts(const int p) { ++(Trace0->pevts.iters[p]); }
inline static long TraceGetIterProcEvts(const int p) { return Trace0->pevts.iters[p]; }
inline static bool TraceRemainsProcEvts(const int p) { return TraceGetIterProcEvts(p)< TraceGetNumProcEvts(p); }

/* pevts-gids */
inline static void TraceSetPtrProcEvtsGids(long **gids) { Trace0->pevts.gids= gids; }
inline static long *TraceGetPtrProcEvtsGids(const int p) { return Trace0->pevts.gids[p]; }
/* inline static long TraceGetProcEvtsGid(const int p, const long it) { return Trace0->pevts.gids[p][it]; } */
/* inline static long TraceGetCurrProcEvtsGid(const int p) { return TraceGetProcEvtsGid(p, TraceGetIterProcEvts(p)); } */
inline static void TraceSetProcEvtGidTo(const int p, const long it,
                                        const long gid) { Trace0->pevts.gids[p][it]= gid; }
inline static void TraceSetCurrProcEvtGidTo(const int p, const long gid) { TraceSetProcEvtGidTo(p, TraceGetIterProcEvts(p), gid); }
inline static void TraceSetCurrProcEvtCurrGid(const int p) { TraceSetCurrProcEvtGidTo(p, TraceGetIterEvts()); }
inline static long TraceGetProcEvtGid(const int p, const long it) { return Trace0->pevts.gids[p][it]; }
inline static long TraceGetGidCurrProcEvt(const int p) { return TraceGetProcEvtGid(p, TraceGetIterProcEvts(p)); }

/* pevts-at */
inline static double TraceGetAtProcEvt(const int p, const long ix) { return TraceGetEvtAt(TraceGetProcEvtGid(p, ix)); }
inline static double TraceGetProcNextEvtDelay(const int p, const long ix) { return TraceGetAtProcEvt(p, ix+ 1)- TraceGetAtProcEvt(p, ix); }
inline static double TraceGetAtCurrProcEvt(const int p) { return TraceGetAtProcEvt(p, TraceGetIterProcEvts(p)); }
inline static double TraceGetAtPrevProcEvt(const int p) { return 0== TraceGetIterProcEvts(p)? TraceGetProcStartTime(p): TraceGetAtProcEvt(p, TraceGetIterProcEvts(p)- 1); }
inline static double TraceGetAtNextProcEvt(const int p) { return (TraceGetNumProcEvts(p)- 1)== TraceGetIterProcEvts(p)? TraceGetProcEndTime(p): TraceGetAtProcEvt(p, TraceGetIterProcEvts(p)+ 1); }
/* pevts-id */
inline static int TraceGetIdProcEvt(const int p, const long it) { return TraceGetEvtId(TraceGetProcEvtGid(p, it)); }
inline static int TraceGetIdCurrProcEvt(const int p) { return TraceGetIdProcEvt(p, TraceGetIterProcEvts(p)); }
inline static int TraceGetIdPrevProcEvt(const int p) { return 0== TraceGetIterProcEvts(p)? -99: TraceGetIdProcEvt(p, TraceGetIterProcEvts(p)- 1); }
inline static int TraceGetIdNextProcEvt(const int p) { return TraceGetIdProcEvt(p, TraceGetIterProcEvts(p)+ 1); }
/* pevts-crit */
inline static void TraceSetCritCurrProcEvt(const int p, const double crit) { TraceSetEvtCrit(TraceGetGidCurrProcEvt(p), crit); }
inline static double TraceGetCritProcEvt(const int p, const long ix) { return TraceGetEvtCrit(TraceGetProcEvtGid(p, ix)); }
inline static double TraceGetProcNextEvtCritDelay(const int p, const long ix) { return TraceGetCritProcEvt(p, ix+ 1)- TraceGetCritProcEvt(p, ix); }
inline static double TraceGetCritCurrProcEvt(const int p) { return TraceGetCritProcEvt(p, TraceGetIterProcEvts(p)); }
inline static double TraceGetCritPrevProcEvt(const int p) { return 0== TraceGetIterProcEvts(p)? TraceGetProcStartTime(p): TraceGetCritProcEvt(p, TraceGetIterProcEvts(p)- 1); }
inline static double TraceGetCritNextProcEvt(const int p) { return (TraceGetNumProcEvts(p)- 1)== TraceGetIterProcEvts(p)? TraceGetCritProcEvt(p, TraceGetNumProcEvts(p)- 1): TraceGetCritProcEvt(p, TraceGetIterProcEvts(p)+ 1); }

/* pevts-slist */
inline static void TraceResetProcEvtSends(const int p, const int i,
                                          const long it) { TraceSetEvtSends(i, TraceGetProcEvtGid(p, it), NULL); }
inline static void TraceResetCurrProcEvtSends(const int p, const int i) { TraceResetProcEvtSends(p, i, TraceGetIterProcEvts(p)); }
inline static IndexList *TraceGetProcEvtSends(const int p, const int i,
                                              const long it) { return TraceGetEvtSends(i, TraceGetProcEvtGid(p, it)); }
inline static IndexList *TraceGetCurrProcEvtSends(const int p, const int i) { return TraceGetProcEvtSends(p, i, TraceGetIterProcEvts(p)); }
inline static bool TraceExistSendsOnCurrProcEvt(const int p, const int i) { return NULL!= TraceGetCurrProcEvtSends(p, i); }
inline static bool TraceNoSendOnCurrProcEvt(const int p) { return !TraceExistSendsOnCurrProcEvt(p, 0)&& !TraceExistSendsOnCurrProcEvt(p, 1); }

/* pevts-rlist */
inline static void TraceResetProcEvtRecvs(const int p, const int i,
                                          const long it) { TraceSetEvtRecvs(i, TraceGetProcEvtGid(p, it), NULL); }
inline static void TraceResetCurrProcEvtRecvs(const int p, const int i) { TraceResetProcEvtRecvs(p, i, TraceGetIterProcEvts(p)); }
inline static IndexList *TraceGetProcEvtRecvs(const int p, const int i,
                                              const long it) { return TraceGetEvtRecvs(i, TraceGetProcEvtGid(p, it)); }
inline static IndexList *TraceGetCurrProcEvtRecvs(const int p, const int i) { return TraceGetProcEvtRecvs(p, i, TraceGetIterProcEvts(p)); }
inline static bool TraceExistRecvsOnCurrProcEvt(const int p, const int i) { return NULL!= TraceGetCurrProcEvtRecvs(p, i); }
inline static bool TraceNoRecvOnCurrProcEvt(const int p) { return !TraceExistRecvsOnCurrProcEvt(p, 0)&& !TraceExistRecvsOnCurrProcEvt(p, 1); }

/* pevts-multi */
inline static bool TraceNoMsgOnCurrProcEvt(const int p) { return TraceNoSendOnCurrProcEvt(p)&& TraceNoRecvOnCurrProcEvt(p); }
/* inline static bool TraceNoMsgOnNextProcEvt(const int p) { return TraceNoSendOnNextProcEvt(p)&& TraceNoRecvOnNextProcEvt(p); } */

/* evts-multi */
inline static void TraceRegisterProcEvt(const int p, const double t,
                                        const int evtId)
{
  TraceSetCurrEvtAt(t);
  TraceSetCurrEvtId(evtId);
  TraceSetCurrEvtProc(p);

  TraceSetCurrProcEvtCurrGid(p);
  TraceIncrIterProcEvts(p);

  TraceIncrIterEvts();
}

/* pevts-time-based-query */
/* assumes increasing order of times at events' log */
#if defined (_EVT_SEARCH_BIN)
inline static long TraceSearchIterProcEvtAt(const int p, const double at,
                                            const long it0, const long it1)
{
  /* this search module is not updated with the nearby event time matching
   * it still follows exact matching - useless for simulated traces
   */
  if(it1- it0< 2) {
    if(fabs(at- TraceGetAtProcEvt(p, it0))< 0.1) {
      return it0;
    } else if(fabs(at- TraceGetAtProcEvt(p, it1))< 0.1) {
      return it1;
    } else {
      printf("ERROR: I shouldn't be here (0)!\n");
      return -1;
    }
  }

  const long itby2= (it0+ it1)/ 2;
  const double timeAtBisect= TraceGetAtProcEvt(p, itby2);
  if(fabs(at- timeAtBisect)< 0.1) {
    return itby2;
  } else if(at< timeAtBisect) {
    return TraceSearchIterProcEvtAt(p, at, it0, itby2- 1);
  } else {
    return TraceSearchIterProcEvtAt(p, at, itby2+ 1, it1);
  }

  /* one should never be here */
  printf("ERROR: I shouldn't be here (1)!\n");
  return -1;
}
#else
inline static long TraceSearchIterProcEvtAt(const int p, const double at,
                                            const long it0, const long it1)
{
  if(it0== it1) {
    return it0;
  }
  const double t0= TraceGetAtProcEvt(p, it0);
  const double t1= TraceGetAtProcEvt(p, it1);
  if(at- t0> t1- at) {
    for(long it= it1; it>= it0; --it) {
      if(TraceGetAtProcEvt(p, it)- at< 0.1) {
        return it;
      }
    }
  } else {
    for(long it= it0; it<= it1; ++it) {
      if(at- TraceGetAtProcEvt(p, it)< 0.1) {
        return it;
      }
    }
  }

  /* one should never be here */
  printf("ERROR: I shouldn't be here (2)!\n");
  return -1;
}
#endif
inline static long TraceGetIterProcEvtAt(const int p, const double at,
                                         const long tiebreakguide)
{
  static int lastp= -1; static double lastat= -1.0; static long lastit= -1;

  long it0= 0, it1= TraceGetNumProcEvts(p)- 1;
  if(p== lastp&& -1!= lastit) {
    const double diff= at- lastat;
    if(fabs(diff)< 0.1) {
      it0= lastit;
      it1= lastit;
    } else if(diff> 0.9) {
      it0= lastit;
    } else if(diff< -0.9) {
      it1= lastit;
    }
  }
  lastp= p; lastat= at;
  long it= TraceSearchIterProcEvtAt(p, at, it0, it1);

  /* 'tiebreakguide': tie-break instantaneous events (usually simulated trace) */
  if(0!= it&& it+ tiebreakguide>= 0&& it+ tiebreakguide< TraceGetNumProcEvts(p)&&
     fabs(at- TraceGetAtProcEvt(p, it+ tiebreakguide))< 0.1) {
    it+= tiebreakguide;
  }
  lastit= it;
  return it;
}
inline static long TraceGetGidProcEvtAt(const int p, const double at,
                                        const long tiebreakguide)
{
  return TraceGetProcEvtGid(p, TraceGetIterProcEvtAt(p, at, tiebreakguide));
}

inline static long TraceSearchIterProcEvtId(const int p, const int id)
{
  for(long ix= TraceGetIterProcEvts(p); ix< TraceGetNumProcEvts(p); ++ix) {
    if(id== TraceGetIdProcEvt(p, ix)) {
      return ix;
    }
  }
  return -1;
}
inline static int TraceGetIdProcEvtAt(const int p, const double at,
                                      const long tiebreakguide)
{
  long it= TraceGetIterProcEvtAt(p, at, tiebreakguide);
  if(-1== it) {
    return -12345;
  }

  if(0== TraceGetIdProcEvt(p, it)) { /* end of MPI call */
    /* return negative to indicate exit of an MPI call (start of useful) */
    if(it> 0) {
      return -TraceGetIdProcEvt(p, it- 1);
    } else {
      return -0;
    }
  }
  return TraceGetIdProcEvt(p, it);
}

/* msgs-num */
inline static void TraceSetNumMsgs(const long num) { Trace0->msgs.num= num; }
inline static long TraceGetNumMsgs() { return Trace0->msgs.num; }

/* msgs-iter */
inline static void TraceResetIterMsgs() { Trace0->msgs.iter= 0; }
inline static void TraceIncrIterMsgs() { ++(Trace0->msgs.iter); }
inline static long TraceGetIterMsgs() { return Trace0->msgs.iter; }

/* msgs-send-at */
inline static void TraceSetPtrMsgsSendAt(double (*at)[3]) { Trace0->msgs.st= at; }
inline static double (*TraceGetPtrMsgsSendAt())[3] { return Trace0->msgs.st; }
inline static void TraceSetMsgSendAt(const long it, const int i,
                                     const double t) { Trace0->msgs.st[it][i]= t; }
inline static double TraceGetMsgSendAt(const long it, const int i) { return Trace0->msgs.st[it][i]; }
inline static double *TraceGetPtrMsgSendAt(const long it, const int i) { return Trace0->msgs.st[it]+ i; }
inline static void TraceSetCurrMsgSendAt(const int i, const double t) { TraceSetMsgSendAt(TraceGetIterMsgs(), i, t); }
inline static double *TraceGetPtrCurrMsgSendAt(const int i) { return TraceGetPtrMsgSendAt(TraceGetIterMsgs(), i); }

/* msgs-recv-at */
inline static void TraceSetPtrMsgsRecvAt(double (*at)[3]) { Trace0->msgs.rt= at; }
inline static double (*TraceGetPtrMsgsRecvAt())[3] { return Trace0->msgs.rt; }
inline static void TraceSetMsgRecvAt(const long it, const int i,
                                     const double t) { Trace0->msgs.rt[it][i]= t; }
inline static double TraceGetMsgRecvAt(const long it, const int i) { return Trace0->msgs.rt[it][i]; }
inline static double *TraceGetPtrMsgRecvAt(const long it, const int i) { return Trace0->msgs.rt[it]+ i; }
inline static void TraceSetCurrMsgRecvAt(const int i, const double t) { TraceSetMsgRecvAt(TraceGetIterMsgs(), i, t); }
inline static double *TraceGetPtrCurrMsgRecvAt(const int i) { return TraceGetPtrMsgRecvAt(TraceGetIterMsgs(), i); }

/* msgs-send-rank */
inline static void TraceSetPtrMsgsSendRank(int *rank) { Trace0->msgs.srank= rank; }
inline static int *TraceGetPtrMsgsSendRank() { return Trace0->msgs.srank; }
inline static int *TraceGetPtrMsgSendRank(const long it) { return Trace0->msgs.srank+ it; }
inline static int *TraceGetPtrCurrMsgSendRank() { return TraceGetPtrMsgSendRank(TraceGetIterMsgs()); }
inline static void TraceSetMsgSendRank(const long it, const int rank) { Trace0->msgs.srank[it]= rank; }
inline static void TraceSetCurrMsgSendRank(const int rank) { TraceSetMsgSendRank(TraceGetIterMsgs(), rank); }
inline static int TraceGetMsgSendRank(const long it) { return Trace0->msgs.srank[it]; }
inline static int TraceGetCurrMsgSendRank() { return TraceGetMsgSendRank(TraceGetIterMsgs()); }

/* msgs-recv-rank */
inline static void TraceSetPtrMsgsRecvRank(int *rank) { Trace0->msgs.rrank= rank; }
inline static int *TraceGetPtrMsgsRecvRank() { return Trace0->msgs.rrank; }
inline static int *TraceGetPtrMsgRecvRank(const long it) { return Trace0->msgs.rrank+ it; }
inline static int *TraceGetPtrCurrMsgRecvRank() { return TraceGetPtrMsgRecvRank(TraceGetIterMsgs()); }
inline static void TraceSetMsgRecvRank(const long it, const int rank) { Trace0->msgs.rrank[it]= rank; }
inline static void TraceSetCurrMsgRecvRank(const int rank) { TraceSetMsgRecvRank(TraceGetIterMsgs(), rank); }
inline static int TraceGetMsgRecvRank(const long it) { return Trace0->msgs.rrank[it]; }
inline static int TraceGetCurrMsgRecvRank() { return TraceGetMsgRecvRank(TraceGetIterMsgs());}

/* msgs-size */
inline static void TraceSetPtrMsgsSize(double *size) { Trace0->msgs.size= size; }
inline static double *TraceGetPtrMsgsSize() { return Trace0->msgs.size; }
inline static double *TraceGetPtrMsgSize(const long it) { return Trace0->msgs.size+ it; }
inline static void TraceSetCurrMsgSize(const double size) { Trace0->msgs.size[TraceGetIterMsgs()]= size; }
inline static double *TraceGetPtrCurrMsgSize() { return TraceGetPtrMsgSize(TraceGetIterMsgs()); }
inline static double TraceGetMsgSize(const long it) { return Trace0->msgs.size[it]; }

/* msgs-tag */
inline static void TraceSetPtrMsgsTag(int *tag) { Trace0->msgs.tag= tag; }
inline static int *TraceGetPtrMsgsTag() { return Trace0->msgs.tag; }
inline static int *TraceGetPtrMsgTag(const long it) { return Trace0->msgs.tag+ it; }
inline static void TraceSetCurrMsgTag(const int tag) { Trace0->msgs.tag[TraceGetIterMsgs()]= tag; }
inline static int *TraceGetPtrCurrMsgTag() { return TraceGetPtrMsgTag(TraceGetIterMsgs()); }

/* msgs-ip[s/r]mem */
inline static IndexList *TraceGetPtrMsgsProcSendIndex() { return Trace0->msgs.ipsmem; }
inline static IndexList *TraceGetPtrMsgsProcRecvIndex() { return Trace0->msgs.iprmem; }

/* psends-nums */
inline static long *TraceGetPtrNumProcSends() { return Trace0->psends.nums; }
inline static void TraceIncrNumProcSends(const int p) { ++(Trace0->psends.nums[p]); }
inline static long TraceGetNumProcSends(const int p) { return Trace0->psends.nums[p]; }

/* psends-iters */
inline static long *TraceGetPtrIterProcSends() { return Trace0->psends.iters; }
inline static void TraceResetItersProcSends() { memset(TraceGetPtrIterProcSends(), 0, sizeof(long)* TraceGetNumProcs()); }
inline static void TraceIncrIterProcSends(const int p) { ++(Trace0->psends.iters[p]); }
inline static long TraceGetIterProcSends(const int p) { return Trace0->psends.iters[p]; }

/* psends-gids */
inline static void TraceSetPtrProcSendsGids(long **gids) { Trace0->psends.gids= gids; }
inline static long *TraceGetPtrProcSendsGids(const int p) { return Trace0->psends.gids[p]; }
inline static void TraceSetProcSendGidTo(const int p, const long it,
                                         const long gid) { Trace0->psends.gids[p][it]= gid; }
inline static void TraceSetCurrProcSendGidTo(const int p, const long gid) { TraceSetProcSendGidTo(p, TraceGetIterProcSends(p), gid); }
inline static long TraceGetProcSendGid(const int p, const long it) { return Trace0->psends.gids[p][it]; }

/* precvs-nums */
inline static long *TraceGetPtrNumProcRecvs() { return Trace0->precvs.nums; }
inline static void TraceIncrNumProcRecvs(const int p) { ++(Trace0->precvs.nums[p]); }
inline static long TraceGetNumProcRecvs(const int p) { return Trace0->precvs.nums[p]; }

/* precvs-iters */
inline static long *TraceGetPtrIterProcRecvs() { return Trace0->precvs.iters; }
inline static void TraceResetItersProcRecvs() { memset(TraceGetPtrIterProcRecvs(), 0, sizeof(long)* TraceGetNumProcs()); }
inline static void TraceIncrIterProcRecvs(const int p) { ++(Trace0->precvs.iters[p]); }
inline static long TraceGetIterProcRecvs(const int p) { return Trace0->precvs.iters[p]; }

/* precvs-gids */
inline static void TraceSetPtrProcRecvsGids(long **gids) { Trace0->precvs.gids= gids; }
inline static long *TraceGetPtrProcRecvsGids(const int p) { return Trace0->precvs.gids[p]; }
inline static void TraceSetProcRecvGidTo(const int p, const long it,
                                         const long gid) { Trace0->precvs.gids[p][it]= gid; }
inline static void TraceSetCurrProcRecvGidTo(const int p, const long gid) { TraceSetProcRecvGidTo(p, TraceGetIterProcRecvs(p), gid); }
inline static long TraceGetProcRecvGid(const int p, const long it) { return Trace0->precvs.gids[p][it]; }

/* psends/precvs-remote */
inline static int TraceGetProcSendRemote(const int p, const long it) { return TraceGetMsgRecvRank(TraceGetProcSendGid(p, it)); }
inline static int TraceGetCurrProcSendRemote(const int p) { return TraceGetProcSendRemote(p, TraceGetIterProcSends(p)); }
inline static int TraceGetProcRecvRemote(const int p, const long it) { return TraceGetMsgSendRank(TraceGetProcRecvGid(p, it)); }
inline static int TraceGetCurrProcRecvRemote(const int p) { return TraceGetProcRecvRemote(p, TraceGetIterProcRecvs(p)); }

/* psend/precv-remote-at */
inline static double *TraceGetProcSendRemoteRecvAts(const int p,
                                                    const long ix) { return TraceGetPtrMsgRecvAt(TraceGetProcSendGid(p, ix), 0); }
inline static double *TraceGetProcRecvRemoteSendAts(const int p,
                                                    const long ix) { return TraceGetPtrMsgSendAt(TraceGetProcRecvGid(p, ix), 0); }

/* psend/precv-at */
inline static void TraceSetProcSendAt(const int p, const long it, const int i,
                                      const double t) { TraceSetMsgSendAt(TraceGetProcSendGid(p, it), i, t); }
inline static double TraceGetProcSendAt(const int p, const long it,
                                        const int i) { return TraceGetMsgSendAt(TraceGetProcSendGid(p, it), i); }
inline static double *TraceGetRecvAtsProcSend(const int p, const long ix) { return TraceGetPtrMsgRecvAt(TraceGetProcSendGid(p, ix), 0); }
inline static double TraceGetRecvAtProcSend(const int p, const long it,
                                            const int i) { return TraceGetMsgRecvAt(TraceGetProcSendGid(p, it), i); }
inline static void TraceSetProcRecvAt(const int p, const long it, const int i,
                                      const double t) { TraceSetMsgRecvAt(TraceGetProcRecvGid(p, it), i, t); }
inline static double TraceGetProcRecvAt(const int p, const long it,
                                        const int i) { return TraceGetMsgRecvAt(TraceGetProcRecvGid(p, it), i); }
inline static double *TraceGetSendAtsProcRecv(const int p, const long ix) { return TraceGetPtrMsgSendAt(TraceGetProcRecvGid(p, ix), 0); }
inline static double TraceGetSendAtProcRecv(const int p, const long it,
                                            const int i) { return TraceGetMsgSendAt(TraceGetProcRecvGid(p, it), i); }
inline static double *TraceGetProcSendAts(const int p, const long it) { return TraceGetPtrMsgSendAt(TraceGetProcSendGid(p, it), 0); }
inline static double *TraceGetProcRecvAts(const int p, const long it) { return TraceGetPtrMsgRecvAt(TraceGetProcRecvGid(p, it), 0); }

/* psend-size */
inline static double TraceGetProcSendSize(const int p, const long it) { return TraceGetMsgSize(TraceGetProcSendGid(p, it)); }

/* psend-precv-multi */
inline static void TraceSetProcSendRecvGids(const int psend, const int precv)
{
  const long gid= TraceGetIterMsgs();
  TraceSetCurrProcSendGidTo(psend, gid);
  TraceIncrIterProcSends(psend);
  TraceSetCurrProcRecvGidTo(precv, gid);
  TraceIncrIterProcRecvs(precv);
}

/* event to send/recv connection */
extern void TraceConnectEvtsToMsgs();

/* inline static bool TraceHasCurrProcEvtSendStart(const int p) { return TraceGetProcEvtsSends(p, TraceGetGidCurrProcEvt(p)); } */

/* pcolls-allnum */
inline static void TraceSetNumAllProcColls(const long allnum) { Trace0->pcolls.allnum= allnum; }
inline static long TraceGetNumAllProcColls() { return Trace0->pcolls.allnum; }

/* pcolls-nums */
inline static long *TraceGetPtrNumProcColls() { return Trace0->pcolls.nums; }
inline static void TraceIncrNumProcColls(const int p) { ++(Trace0->pcolls.nums[p]); }
inline static long TraceGetNumProcColls(const int p) { return Trace0->pcolls.nums[p]; }

/* pcolls-iters */
inline static long *TraceGetPtrIterProcColls() { return Trace0->pcolls.iters; }
inline static void TraceResetItersProcColls() { memset(TraceGetPtrIterProcColls(), 0, sizeof(long)* TraceGetNumProcs()); }
inline static void TraceIncrIterProcColls(const int p) { ++(Trace0->pcolls.iters[p]); }
inline static long TraceGetIterProcColls(const int p) { return Trace0->pcolls.iters[p]; }

/* pcolls-at */
inline static void TraceSetPtrProcCollsAt(double (**at)[2]) { Trace0->pcolls.at= at; }
inline static double (*TraceGetPtrProcCollsAt(const int p))[2] { return Trace0->pcolls.at[p]; }
inline static void TraceSetProcCollAt(const int p, const long it, const int i,
                                      const double at) { Trace0->pcolls.at[p][it][i]= at; }
inline static void TraceSetCurrProcCollAt(const int p, const int i,
                                          const double at) { TraceSetProcCollAt(p, TraceGetIterProcColls(p), i, at); }
inline static double TraceGetAtProcColl(const int p, const long ix,
                                        const int i) { return Trace0->pcolls.at[p][ix][i]; }
inline static double TraceGetAtCurrProcColl(const int p, const int i) { return TraceGetAtProcColl(p, TraceGetIterProcColls(p), i); }

/* pcolls-comm */
inline static void TraceSetPtrProcCollsComm(int **comm) { Trace0->pcolls.comm= comm; }
inline static int *TraceGetPtrProcCollsComm(const int p) { return Trace0->pcolls.comm[p]; }
inline static void TraceSetProcCollComm(const int p, const long ix,
                                        const int comm) { Trace0->pcolls.comm[p][ix]= comm; }
inline static void TraceSetCurrProcCollComm(const int p, const int comm) { TraceSetProcCollComm(p, TraceGetIterProcColls(p), comm); }
inline static int TraceGetCommProcColl(const int p, const long ix) { return Trace0->pcolls.comm[p][ix]; }
inline static int TraceGetCommCurrProcColl(const int p) { return TraceGetCommProcColl(p, TraceGetIterProcColls(p)); }

/* pcolls-multi */
inline static void TraceStartProcColl(const int p, const double t) { TraceSetCurrProcCollAt(p, 0, t); }
inline static void TraceEndProcColl(const int p, const double t) { TraceSetCurrProcCollAt(p, 1, t); TraceIncrIterProcColls(p); }

/* resets */
inline static void TraceResetProcIters()
{
  TraceResetItersProcEvts();
  TraceResetItersProcSends();
  TraceResetItersProcRecvs();
  TraceResetItersProcColls();
}

/* allocation etc. */
extern void TraceAllocAndInitLevel1Data();

#endif  /* REPLAY_TRACE_DATA_H__ */
