/*
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#include"common.h"
#include"utils.h"
#include"paraver_file.h"
#include"trace_data.h"
#include"paraver.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<float.h>
#include<math.h>

/* state/event-identifiers:
 *  0: computation aka useful
 * >0: MPI
 * -1: tracing ended
 * -2: tracing disabled - ideal includes, but not useful
 * -3: flush
 * -4: trace-init
 * -99: invalid
 */

inline static bool isUseful(const int evt) { return 0== evt; }
inline static bool isMPI(const int evt) { return 0< evt; }
inline static bool isDisabledEvt(const int evt) { return -2== evt; }
inline static bool isFlush(const int evt) { return -3== evt; }
inline static bool isTraceInit(const int evt) { return -4== evt; }
inline static bool isUndefined(const int evt) { return -99== evt; }

static struct {
  double *evtAt;
  int *evtId;
  double *tickAt;
  bool *tracingOn;
} last= { NULL, NULL, NULL, NULL };

inline static void initLasts(const int np)
{
  for(int ip= 0; ip< np; ++ip) {
    last.evtAt[ip]= -1.0;
    last.evtId[ip]= -99;
    last.tickAt[ip]= -1.0;
    last.tracingOn[ip]= true; /* appears in trace only when tracing is disabled */
  }
}
inline static void allocLasts(const int np)
{
  last.evtAt= (double *) malloc(sizeof(double)* np);
  last.evtId= (int *) malloc(sizeof(int)* np);
  last.tickAt= (double *) malloc(sizeof(double)* np);
  last.tracingOn= (bool *) malloc(sizeof(bool)* np);
  initLasts(np);
}
inline static void freeLasts()
{
  FREE_IF(last.tracingOn);
  FREE_IF(last.tickAt);
  FREE_IF(last.evtId);
  FREE_IF(last.evtAt);
}
inline static void setLastEvt(const int p, const double t, const int e) { last.evtAt[p]= t; last.evtId[p]= e; }
inline static void setTracingState(const int p, const bool enable) { last.tracingOn[p]= enable; }
inline static bool isTracingDisabled(const int p) { return !last.tracingOn[p]; }

/* count */
inline static void startUseful(const int p, const double at, const int evt)
{
  ErrorIf(isUseful(last.evtId[p])&& fabs(last.evtAt[p]- at)> 0.1,
          "%d: overwrite useful-start %.0lf -> %.0lf\n", p, last.evtAt[p], at);
  setLastEvt(p, at, evt);
}
inline static double usefulBurstLen(const int p, const double at)
{
  if(unlikely(isUndefined(last.evtId[p]))) {
    return 0.0;
  }
  ErrorIf(!isUseful(last.evtId[p])&& fabs(last.evtAt[p]- at)> 0.1,
          "%d: useful-end at %.0lf, last-event %d\n", p, at, last.evtId[p]);

  return at- last.evtAt[p];
}
inline static void startMPI(const int p, const double at, const int evt)
{
  ErrorIf(isMPI(last.evtId[p])&& fabs(last.evtAt[p]- at)> 0.1,
          "%d: overwrite MPI-start %.0lf -> %.0lf\n", p, last.evtAt[p], at);
  setLastEvt(p, at, evt);
}
inline static double mpiBurstLen(const int p, const double at)
{
  if(unlikely(isUndefined(last.evtId[p]))) {
    return 0.0;
  }
  ErrorIf(!isMPI(last.evtId[p])&& fabs(last.evtAt[p]- at)> 0.1,
          "%d: MPI-end at %.0lf, last-event %d\n", p, at, last.evtId[p]);
  return at- last.evtAt[p];
}

static void countMPIEvt(const int p, const double at, const int evt)
{
  switch(evt) {
  case 0:                       /* end MPI, start comp */
    TraceAddProcMPI(p, mpiBurstLen(p, at));
    startUseful(p, at, evt);
    break;
  case 31:
    TraceSetMPIInitEvt(p);
  /* fall through */
  default:                      /* end comp, start MPI */
    if(!isUseful(last.evtId[p])) {
      Debug1("%d: invent useful burst %.0lf -> %.0lf\n", p,
             TraceGetProcStartTime(p), at);
      startUseful(p, TraceGetProcStartTime(p), 0);
    }
    TraceAddProcComp(p, usefulBurstLen(p, at));
    startMPI(p, at, evt);
  }

  TraceIncrNumProcEvts(p);
}

static void countMPICollEvtComm(const int p, const int comm)
{
  if(0== comm) {
    return;  /* most likely at the end */
  }

  if(!ParaverEvtIsColl(last.evtId[p])) {
    ErrorIf(!ParaverEvtIsColl(last.evtId[p]),
            "%d: comm-id without ongoing collective (ongoing: %d-MPI_%s)\n",
            p, last.evtId[p], GetParaverMPIEvtName(last.evtId[p]));
  }

  if(ParaverCollEvtIsDimemasCompliant(last.evtId[p])) {
    TraceIncrNumProcColls(p);
  }
}
static void countAppEvt(const int p, const double at, const int start)
{
  switch(start) {
  case 1:
    if(TraceGetProcStartTime(p)> -0.1&& fabs(TraceGetProcStartTime(p)- at)> 0.1) {
      Debug1("%d: overwrite proc-start (app-event) %.0lf -> %.0lf\n", p,
             TraceGetProcStartTime(p), at);
    }
    TraceStartProcTimeline(p, at);
    startUseful(p, at, 0);
    break;
  case 0:
    if(TraceGetProcEndTime(p)> -0.1&& fabs(TraceGetProcEndTime(p)- at)> 0.1) {
      Debug1("%d: overwrite proc-end (app-event) %.0lf -> %.0lf\n", p,
             TraceGetProcEndTime(p), at);
    }
    TraceEndProcTimeline(p, at);
    if(isUseful(last.evtId[p])) {
      Debug1("%d: end useful (app-event)\n", p);
      TraceAddProcComp(p, usefulBurstLen(p, at));
    }
    setLastEvt(p, at, -1);
    break;
  default:
    Error("Invalid event-value with event-type application (%d)\n", start);
    break;
  }
}
static void countTraceInitEvt(const int p, const double at, const int start)
{
  switch(start) {
  case 1:
    ErrorIf(isTraceInit(last.evtId[p]),
            "%d: start trace-init at %.0lf overwrites ongoing trace-init since %.0lf\n",
            p, at, last.evtAt[p]);
    if(isUseful(last.evtId[p])) {
      TraceAddProcComp(p, usefulBurstLen(p, at));
    }
    setLastEvt(p, at, -4);
    TraceIncrNumProcEvts(p);
    break;
  case 0:
    ErrorIf(!isTraceInit(last.evtId[p]),
            "%d: end trace-init at %.0lf, last event %d at %.0lf\n",
            p, at, last.evtId[p], last.evtAt[p]);
    startUseful(p, at, 0);
    TraceIncrNumProcEvts(p);
    break;
  default:
    Error("%d: invalid trace-init-event value (%d); ignored\n", p, start);
    break;
  }
}
static void countFlushEvt(const int p, const double at, const int start)
{
  switch(start) {
  case 1:
    ErrorIf(isFlush(last.evtId[p]),
            "%d: start flush-event at %.0lf overwrite ongoing flush since %.0lf\n",
            p, at, last.evtAt[p]);
    if(isUseful(last.evtId[p])) {
      TraceAddProcComp(p, usefulBurstLen(p, at));
    }
    TraceStartProcFlush(p, at);
    setLastEvt(p, at, -3);
    TraceIncrNumProcEvts(p);
    break;
  case 0:
    ErrorIf(!isFlush(last.evtId[p]),
            "%d: end flush-event at %.0lf, last event %d at %.0lf\n",
            p, at, last.evtId[p], last.evtAt[p]);
    TraceEndProcFlush(p, at);
    if(isTracingDisabled(p)) {
      setLastEvt(p, at, -2);
    } else {
      startUseful(p, at, 0);
    }
    TraceIncrNumProcEvts(p);
    break;
  default:
    Error("%d: invalid flush-event value (%d); ignore\n", p, start);
    break;
  }
}
static void countTraceability(const int p, const double at, const int start)
{
  switch(start) {
  case 1:
    ErrorIf(!isDisabledEvt(last.evtId[p]),
            "%d: enable tracing at %.0lf overwrites already enabled since %.0lf\n",
            p, last.evtAt[p], -TraceProcDisabledAt(p));
    TraceEndProcDisabled(p, at);
    setTracingState(p, true);
    startUseful(p, at, 0);
    TraceIncrNumProcEvts(p);
    break;
  case 0:
    ErrorIf(isDisabledEvt(last.evtId[p]),
            "%d: disable tracing at %.0lf overwrites already disabled since %.0lf\n",
            p, at, TraceProcDisabledAt(p));
    if(isUseful(last.evtId[p])) {
      TraceAddProcComp(p, usefulBurstLen(p, at));
    }
    TraceStartProcDisabled(p, at);
    setTracingState(p, false);
    setLastEvt(p, at, -2);
    TraceIncrNumProcEvts(p);
    break;
  default:
    Error("%d: invalid traceability-event value (%d); ignore\n", p, start);
    break;
  }
}
static void countEvt(char *const line)
{
  char *ptr= ParaverRecordNextNumNth(line, 3);
  const int p= atoi(ptr)- 1;
  ptr= ParaverRecordNextNumNth(ptr, 2);
  last.tickAt[p]= atof(ptr);

  /* TODO: maybe move this check to a separate reading and finish as soon as everyone is started */
  if(unlikely(TraceGetProcStartTime(p)< 0.0)) {
    Debug1("%d: Starting from trace event at %.0lf ns\n", p, last.tickAt[p]);
    TraceStartProcTimeline(p, last.tickAt[p]);
    startUseful(p, last.tickAt[p], 0);
  }

  ptr= strchr(ptr, ':');
  while(NULL!= ptr) {
    ++ptr;
    long long type= atoll(ptr);
    ptr= ParaverRecordNextNum(ptr);
    switch(type) {
    case 50000001:  /* mpi p2p */
    case 50000002:  /* mpi collective */
    case 50000003:  /* mpi other */
    case 50000004:  /* mpi rma */
    case 50000005:  /* mpi i/o */
      countMPIEvt(p, last.tickAt[p], atoi(ptr));  /* increases #evts */
      break;
    case 50100004:  /* communicator of mpi collective */
      countMPICollEvtComm(p, atoi(ptr));
      break;
    case 40000001:  /* application */
      countAppEvt(p, last.tickAt[p], atoi(ptr));
      break;
    case 40000002:  /* trace init */
      countTraceInitEvt(p, last.tickAt[p], atoi(ptr));
      break;
    case 40000003:  /* flush */
      countFlushEvt(p, last.tickAt[p], atoi(ptr));
      break;
    case 40000012:  /* tracing enabled/disabled */
      countTraceability(p, last.tickAt[p], atoi(ptr));
      break;
    default:
      break;
    }
    ptr= strchr(ptr, ':');
  }
}
static void countMsg(char *const line)
{
  char *ptr= ParaverRecordNextNumNth(line, 3);
  TraceIncrNumProcSends(atoi(ptr)- 1);

  ptr= ParaverRecordNextNumNth(ptr, 6);
  TraceIncrNumProcRecvs(atoi(ptr)- 1);
}
static void evtsAndCommsCounter(char *const line)
{
  if(0== strlen(line)) {
    return;
  }
  switch(atoi(line)) {
  case 2:  /* events */
    countEvt(line);
    break;
  case 3:  /* msgs */
    countMsg(line);
    break;
  default:
    break;
  }
}
static void finishCount(const int np)
{
  double tEnd= 0.0;
  for(int ip= 0; ip< np; ++ip) {
    if(isUseful(last.evtId[ip])) {
      Debug1("%d: end useful (end-of-trace) at %.0lf\n", ip, last.evtAt[ip]);
      TraceAddProcComp(ip, usefulBurstLen(ip, last.tickAt[ip]));
    }

    if(isTracingDisabled(ip)) {
      Debug1("%d: ending tracing-disabled region (end-of-trace) at %.0lf\n",
             ip, TraceProcDisabledAt(ip));
      TraceEndProcDisabled(ip, last.tickAt[ip]);
    }

    tEnd= MAX(MAX(last.tickAt[ip], TraceGetProcEndTime(ip)), tEnd);
  }

  for(int ip= 0; ip< np; ++ip) {
    if(TraceGetProcEndTime(ip)< 0.0) {
      TraceEndProcTimeline(ip, last.tickAt[ip]);
      TraceAddProcComp(ip, tEnd- TraceGetProcEndTime(ip));
    }
  }

  TraceCalculateTimelineExtremes();
}

/* read */
static void readMPIEvt(const int p, const double at, const int evt)
{
  setLastEvt(p, at, evt);
  TraceRegisterProcEvt(p, at, evt);
}
static void readMPICollEvt(const int p, const double t, const int collEvt)
{
  switch(collEvt) {
  case 0:                       /* end of a collective */
    if(likely(ParaverCollEvtIsDimemasCompliant(last.evtId[p]))) {
      TraceEndProcColl(p, t);
    }
    break;
  default:                      /* start of a collective */
    if(likely(ParaverCollEvtIsDimemasCompliant(collEvt))) {
      TraceStartProcColl(p, t);
    }
    break;
  }
}
static void readMPICollEvtComm(const int p, const int comm)
{
  if(0== comm) {
    return;  /* most likely at the end */
  }
  ErrorIf(!ParaverEvtIsColl(last.evtId[p]),
          "%d: comm-id without ongoing collective (ongoing: %d-MPI_%s)\n",
          p, last.evtId[p], GetParaverMPIEvtName(last.evtId[p]));

  if(likely(ParaverCollEvtIsDimemasCompliant(last.evtId[p]))) {
    TraceSetCurrProcCollComm(p, comm);
  }
}
static void readTraceInitEvt(const int p, const double at, const int evt)
{
  switch(evt) {
  case 1:
    Debug1("%d: start trace-init at %.0lf\n", p, at);
    setLastEvt(p, at, -4);
    TraceRegisterProcEvt(p, at, -4);
    break;
  case 0:
    Debug1("%d: end trace-init at %.0lf\n", p, at);
    setLastEvt(p, at, 0);
    TraceRegisterProcEvt(p, at, 0);
    break;
  default:
    Error("%s: unknown trace-init value (%d)\n", __func__, evt);
    break;
  }
}
static void readFlushEvt(const int p, const double at, const int evt)
{
  switch(evt) {
  case 1:
    Debug1("%d: start flush at %.0lf\n", p, at);
    setLastEvt(p, at, -3);
    TraceRegisterProcEvt(p, at, -3);
    break;
  case 0:
    Debug1("%d: end flush at %.0lf\n", p, at);
    if(isTracingDisabled(p)) {
      setLastEvt(p, at, -2);
      TraceRegisterProcEvt(p, at, -2);
    } else {
      setLastEvt(p, at, 0);
      TraceRegisterProcEvt(p, at, 0);
    }
    break;
  default:
    ErrorNL("%s: Unknown flush-event value (%d)\n", __func__, evt);
    break;
  }
}
static void readTraceability(const int p, const double at, const int evt)
{
  switch(evt) {
  case 1:
    Debug1("%d: enable tracing at %.0lf\n", p, at);
    setTracingState(p, true);
    setLastEvt(p, at, 0);
    TraceRegisterProcEvt(p, at, 0);
    break;
  case 0:
    Debug1("%d: tracing disabled at %.0lf\n", p, at);
    setTracingState(p, false);
    setLastEvt(p, at, -2);
    TraceRegisterProcEvt(p, at, -2);
    break;
  default:
    ErrorNL("%s: Unknown traceability-event value (%d)\n", __func__, evt);
    break;
  }
}

static void readEvt(char *const line)
{
  char *ptr= ParaverRecordNextNumNth(line, 3);
  const int p= atoi(ptr)- 1;
  ptr= ParaverRecordNextNumNth(ptr, 2);
  last.tickAt[p]= atof(ptr);
  ptr= strchr(ptr, ':');
  while(NULL!= ptr) {
    ++ptr;
    long long type= atoll(ptr);
    ptr= ParaverRecordNextNum(ptr);
    switch(type) {
    case 50000002:  /* mpi collective */
      readMPICollEvt(p, last.tickAt[p], atoi(ptr));
    /* fall through!! */
    case 50000001:  /* mpi p2p */
    case 50000003:  /* mpi other */
    case 50000004:  /* mpi rma */
    case 50000005:  /* mpi i/o */
      readMPIEvt(p, last.tickAt[p], atoi(ptr));
      break;
    case 50100004:  /* communicator of mpi collective */
      readMPICollEvtComm(p, atoi(ptr)- 1);
      break;
    case 40000001:  /* application */
      break;
    case 40000002:  /* trace init */
      readTraceInitEvt(p, last.tickAt[p], atoi(ptr));
      break;
    case 40000003:  /* flush */
      readFlushEvt(p, last.tickAt[p], atoi(ptr));
      break;
    case 40000012:  /* tracing enabled/disabled */
      readTraceability(p, last.tickAt[p], atoi(ptr));
      break;
    default:
      break;
    }
    ptr= strchr(ptr, ':');
  }
}

static void readMsg(char *const line)
{
  /* 0:1   :2    :3    :4      :5    :6    :7   :8    :9    :10     :11   :12   :13  :14  */
  /* s s    s     x     s       x     x     s    s     x     s       x     x     x    x   */
  /* 3:scpu:stask:srank:sthread:lsend:psend:rcpu:rtask:rrank:rthread:lrecv:precv:size:tag */
  /*                     0   1   2  3   4   5   6   7   8  9   0   1   2   3  4 */
  char *ptr= ParaverRecordNextNumNth(line, 3);
  TraceSetCurrMsgSendRank(atoi(ptr)- 1);

  ptr= ParaverRecordNextNumNth(ptr, 2);
  TraceSetCurrMsgSendAt(0, atof(ptr));

  ptr= ParaverRecordNextNum(ptr);
  TraceSetCurrMsgSendAt(1, atof(ptr));

  ptr= ParaverRecordNextNumNth(ptr, 3);
  TraceSetCurrMsgRecvRank(atoi(ptr)- 1);

  ptr= ParaverRecordNextNumNth(ptr, 2);
  TraceSetCurrMsgRecvAt(0, atof(ptr));

  ptr= ParaverRecordNextNum(ptr);
  TraceSetCurrMsgRecvAt(1, atof(ptr));

  ptr= ParaverRecordNextNum(ptr);
  TraceSetCurrMsgSize(atof(ptr));

  ptr= ParaverRecordNextNum(ptr);
  TraceSetCurrMsgTag(atoi(ptr));

  TraceSetProcSendRecvGids(TraceGetCurrMsgSendRank(),
                           TraceGetCurrMsgRecvRank());

  TraceIncrIterMsgs();
}
static void MPIEvtsAndCommsReader(char *const line)
{
  if(0== strlen(line)) {
    return;
  }
  switch(atoi(line)) {
  case 2:  /* events */
    readEvt(line);
    break;
  case 3:  /* msgs */
    readMsg(line);
    break;
  default:
    break;
  }
}
static void finishRead(const int np)
{
  if(TraceGetNumEvts()!= TraceGetIterEvts()) {
    Error("#events= %ld, event-iter= %ld\n",
          TraceGetNumEvts(), TraceGetIterEvts());
    for(int ip= 0; ip< np; ++ip) {
      if(TraceGetNumProcEvts(ip)!= TraceGetIterProcEvts(ip)) {
        Error("%d: #proc-events= %ld, #proc-event-iter= %ld\n", ip,
              TraceGetNumProcEvts(ip), TraceGetIterProcEvts(ip));
      }
    }
  } else {
    Debug1("event-count check successful\n");
  }

  if(TraceGetNumMsgs()!= TraceGetIterMsgs()) {
    Error("#messages= %ld, message-iter= %ld\n",
          TraceGetNumMsgs(), TraceGetIterMsgs());
    for(int ip= 0; ip< np; ++ip) {
      if(TraceGetNumProcSends(ip)!= TraceGetIterProcSends(ip)) {
        Error("%d: #proc-sends= %ld, #proc-send-event-iter= %ld\n", ip,
              TraceGetNumProcSends(ip), TraceGetIterProcSends(ip));
      }
      if(TraceGetNumProcRecvs(ip)!= TraceGetIterProcRecvs(ip)) {
        Error("%d: #proc-sends= %ld, #proc-send-event-iter= %ld\n", ip,
              TraceGetNumProcRecvs(ip), TraceGetIterProcRecvs(ip));
      }
    }
  } else {
    Debug1("message-count check successful\n");
  }

  bool mismatch= false;
  for(int ip= 0; ip< np; ++ip) {
    if(TraceGetNumProcColls(ip)!= TraceGetIterProcColls(ip)) {
      mismatch= true;
      Error("%d: #proc-collectives= %ld, #proc-collective-iter= %ld\n", ip,
            TraceGetNumProcColls(ip), TraceGetIterProcColls(ip));
    }
  }
  if(!mismatch) {
    Debug1("collective-count check successful\n");
  }

  if(GlOpts.evt_mon.enabled) {
    printf("Event-driven monitoring rank:");
    if(GlOpts.evt_mon.rank< 0|| GlOpts.evt_mon.rank>= TraceGetNumProcs()) {
      int pMaxUseful= 0;
      for(int ip= 1; ip< TraceGetNumProcs(); ++ip) {
        if(TraceGetProcCompDuration(ip)- TraceGetProcCompDuration(pMaxUseful)> 0.1) {
          pMaxUseful= ip;
        }
      }

      printf(" %d ->", GlOpts.evt_mon.rank);
      GlOpts.evt_mon.rank= pMaxUseful;
    }
    printf(" %d\n", GlOpts.evt_mon.rank);
  }
}

inline static int determineHeaderLens(const int n, const double *const avg,
                                      const char **names, int *lens, const int np)
{
  for(int i= 0; i< n; ++i) {
    lens[i]= strlen(names[i]);
    if(avg[i]> 9.0) {
      const int x= ((int) floor(log10(avg[i])))+ 4;
      lens[i]= MAX(x, lens[i]);
    }
  }
  const int ranklen= (int) floor(log10((double) np))+ 1;
  return MAX(ranklen, 4);
}
static void showAggregated(const int np)
{
  if(!GlOpts.show_opts.profile) {
    return;
  }
#define NAGG 11
  double avg[NAGG];
  memset(avg, 0, sizeof(double)* NAGG);

  const double n2us= 1.0e-3;
  for(int ip= 0; ip< np; ++ip) {
    avg[0]+= (((double) TraceGetRuntime())- TraceGetProcMPIDuration(ip))* n2us;
    avg[1]+= TraceGetProcMPIDuration(ip)* n2us;
    avg[2]+= TraceGetProcCompDuration(ip)* n2us;
    avg[3]+= TraceGetProcFlushDuration(ip)* n2us;
    avg[4]+= TraceGetProcDisabledDuration(ip)* n2us;
    avg[5]+= TraceGetProcStartTime(ip)* n2us;
    avg[6]+= TraceGetProcEndTime(ip)* n2us;
    avg[7]+= (double) TraceGetNumProcEvts(ip);
    avg[8]+= (double) TraceGetNumProcSends(ip);
    avg[9]+= (double) TraceGetNumProcRecvs(ip);
    avg[10]+= (double) TraceGetNumProcColls(ip);
  }

  FILE *fp= NULL;
  {
    const int len= strlen(GlOpts.filename)+ 26;
    char *fn= (char *) malloc(sizeof(char)* len);
    memset(fn, 0, sizeof(char)* len);
    memcpy(fn, GlOpts.filename, len- 26);
    /* printf("fn0= \"%s\"\n", fn); */
    strcat(fn, ".clocktalk.aggregated.txt");
    /* printf("fn1= \"%s\", fn[-3]= '%c', fn[-2]= '%c', fn[-1]= '%c'\n", fn, fn[len- 3], fn[len- 2], fn[len- 1]); */
    fp= fopen(fn, "w");
    FREE_IF(fn);
  }
  const char *hnames[NAGG]= { "Outside-MPI", "MPI", "Useful", "Flush", "Disabled", "start", "end", "nevts", "nsends", "nrecvs", "ncolls" };
  int hlens[NAGG]= { 0 };
  const int ranklen= determineHeaderLens(NAGG, avg, hnames, hlens, np);
  fprintf(fp, "%*s", ranklen, "rank");
  for(int i= 0; i< NAGG; ++i) {
    fprintf(fp, " %*s", hlens[i], hnames[i]);
  }
  fprintf(fp, "\n");

  double min[NAGG], max[NAGG];
  for(int i= 0; i< NAGG; ++i) {
    min[i]= DBL_MAX;
  }
  memset(max, 0, sizeof(double)* NAGG);

  for(int ip= 0; ip< np; ++ip) {
    const double oom= (((double) TraceGetRuntime())- TraceGetProcMPIDuration(
                         ip))* n2us;
    min[0]= MIN(oom, min[0]); max[0]= MAX(oom, max[0]);

    const double mpi= TraceGetProcMPIDuration(ip)* n2us;
    min[1]= MIN(mpi, min[1]); max[1]= MAX(mpi, max[1]);

    const double use= TraceGetProcCompDuration(ip)* n2us;
    min[2]= MIN(use, min[2]); max[2]= MAX(use, max[2]);

    const double flush= TraceGetProcFlushDuration(ip)* n2us;
    min[3]= MIN(flush, min[3]); max[3]= MAX(flush, max[3]);

    /* const double disable= 0.0; */
    const double disable= TraceGetProcDisabledDuration(ip)* n2us;
    min[4]= MIN(disable, min[4]); max[4]= MAX(disable, max[4]);

    const double t0= TraceGetProcStartTime(ip)* n2us;
    min[5]= MIN(t0, min[5]); max[5]= MAX(t0, max[5]);

    const double t1= TraceGetProcEndTime(ip)* n2us;
    min[6]= MIN(t1, min[6]); max[6]= MAX(t1, max[6]);

    const double nevts= (double) TraceGetNumProcEvts(ip);
    min[7]= MIN(nevts, min[7]); max[7]= MAX(nevts, max[7]);

    const double nsends= (double) TraceGetNumProcSends(ip);
    min[8]= MIN(nsends, min[8]); max[8]= MAX(nsends, max[8]);

    const double nrecvs= (double) TraceGetNumProcRecvs(ip);
    min[9]= MIN(nrecvs, min[9]); max[9]= MAX(nrecvs, max[9]);

    const double ncolls= (double) TraceGetNumProcColls(ip);
    min[10]= MIN(ncolls, min[10]); max[10]= MAX(ncolls, max[10]);

    fprintf(fp,
            "%*d %*.2lf %*.2lf %*.2lf %*.2lf %*.2lf %*.2lf %*.2lf %*.0lf %*.0lf %*.0lf %*.0lf\n",
            ranklen, ip, hlens[0], oom, hlens[1], mpi, hlens[2], use, hlens[3], flush,
            hlens[4], disable, hlens[5], t0, hlens[6], t1, hlens[7], nevts,
            hlens[8], nsends, hlens[9], nrecvs, hlens[10], ncolls);
  }

  fprintf(fp, "\n%*s", ranklen, "");
  for(int i= 0; i< NAGG; ++i) {
    fprintf(fp, " %*s", hlens[i], hnames[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "%*s", ranklen, "Tot");
  for(int i= 0; i< NAGG; ++i) {
    fprintf(fp, " %*.2lf", hlens[i], avg[i]);
  }
  fprintf(fp, "\n");

  const double pfactor= 1.0/ ((double) np);
  for(int i= 0; i< NAGG; ++i) {
    avg[i]*= pfactor;
  }

  fprintf(fp, "%*s", ranklen, "Avg");
  for(int i= 0; i< NAGG; ++i) {
    fprintf(fp, " %*.2lf", hlens[i], avg[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "%*s", ranklen, "Max");
  for(int i= 0; i< NAGG; ++i) {
    fprintf(fp, " %*.2lf", hlens[i], max[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "%*s", ranklen, "Min");
  for(int i= 0; i< NAGG; ++i) {
    fprintf(fp, " %*.2lf", hlens[i], min[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "%*s", ranklen, "LB");
  for(int i= 0; i< NAGG; ++i) {
    fprintf(fp, " %*.0lf%%", hlens[i]- 1, avg[i]* 100.0/ max[i]);
  }
  fprintf(fp, "\n");
#undef NAGG
}

inline static double processParaverFile(ParaverFile *const file,
                                        void (*processor)(char *const))
{
  ParaverFileSetLineProcessor(file, processor);
  ParaverFileReloadRecords(file);
  return ParaverFileProcess(file, !GlOpts.show_opts.timings);
}
int ReadParaverFile(const char *const fn)
{
  ParaverFile *file= ParaverFileOpen(fn);

  SetWorkingTrace(CreateTrace(file));

  const int np= ParaverFileGetNumProcs(file);
  allocLasts(np);

  const double tioCount= processParaverFile(file, evtsAndCommsCounter);
  finishCount(np);              /* closes states */

  TraceAllocAndInitLevel1Data();

  TraceResetIterEvts();
  TraceResetIterMsgs();
  TraceResetItersProcEvts();
  TraceResetItersProcSends();
  TraceResetItersProcRecvs();
  TraceResetItersProcColls();

  initLasts(np);

  const double tioRead= processParaverFile(file, MPIEvtsAndCommsReader);
  finishRead(np);               /* match iter and nums */

  showAggregated(np);

  freeLasts();

  ParaverFileClose(file); file= NULL;
  Debug1("File I/O: during count= %.1lfs, during read= %.1lfs\n", tioCount,
         tioRead);

  return 0;
}

const char *GetParaverMPIEvtName(const int ev)
{
  switch(ev) {
  case -99:
    return "Invalid";
  case -4:
    return "Trace-init";
  case -3:
    return "Flush";
  case -2:
    return "Disable-tracing";
  case -1:
    return "End-tracing";
  default:
    return ParaverFileGetMPIName(ev);
  }

  return "Strange";
}
