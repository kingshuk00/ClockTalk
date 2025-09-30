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
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>

inline static void allocComms(TraceData *const t, const int numComms,
                              const long sizeAllComms)
{
  t->comms.num= numComms;
  t->comms.sizes= (int *) malloc(sizeof(int)* numComms);
  memset(t->comms.sizes, 0, sizeof(int)* numComms);
  t->comms.ranks= (int **) malloc(sizeof(int *)* numComms);
  t->comms.ranks[0]= (int *) malloc(sizeof(int)* sizeAllComms);
  memset(t->comms.ranks[0], 0, sizeof(int)* sizeAllComms);
}
inline static void allocLevel0Data(TraceData *const t)
{
  const int np= t->numprocs;
  t->timeline.extents= (double (*)[2]) malloc(sizeof(double[2])* np);
  t->timeline.tcomp= (double *) malloc(sizeof(double)* np);
  t->timeline.tmpi= (double *) malloc(sizeof(double)* np);
  t->timeline.tflush= (double *) malloc(sizeof(double)* np);
  t->timeline.tdisabled= (double *) malloc(sizeof(double)* np);
  t->timeline.disabledAt= (double *) malloc(sizeof(double)* np);
  t->timeline.hasTraceInit= (bool *) malloc(sizeof(bool)* np);
  t->timeline.hasMPIInit= (bool *) malloc(sizeof(bool)* np);
  const size_t size= sizeof(long)* np;
  t->pevts.nums= (long *) malloc(size);
  t->pevts.iters= (long *) malloc(size);
  t->pevts.gids= (long **) malloc(sizeof(long *)* np);
  t->psends.nums= (long *) malloc(size);
  t->psends.iters= (long *) malloc(size);
  t->psends.gids= (long **) malloc(sizeof(long *)* np);
  t->precvs.nums= (long *) malloc(size);
  t->precvs.iters= (long *) malloc(size);
  t->precvs.gids= (long **) malloc(sizeof(long *)* np);
  t->pcolls.nums= (long *) malloc(size);
  t->pcolls.iters= (long *) malloc(size);
  t->pcolls.at= (double (**)[2]) malloc(sizeof(double (*)[2])* np);
  t->pcolls.comm= (int **) malloc(sizeof(int *)* np);
}
inline static void initLevel0Data(TraceData *const t)
{
  const int np= t->numprocs;
  for(int ip= 0; ip< np; ++ip) {
    t->timeline.extents[ip][0]= -1.0;
    t->timeline.extents[ip][1]= -1.0;
    t->timeline.disabledAt[ip]= -1.0;
  }
  memset(t->timeline.tcomp, 0, sizeof(double)* np);
  memset(t->timeline.tmpi, 0, sizeof(double)* np);
  memset(t->timeline.tflush, 0, sizeof(double)* np);
  memset(t->timeline.tdisabled, 0, sizeof(double)* np);
  /* memset(t->timeline.disabledAt, 0, sizeof(double)* np); */
  memset(t->timeline.hasTraceInit, 0, sizeof(bool)* np);
  memset(t->timeline.hasMPIInit, 0, sizeof(bool)* np);
  const size_t size= sizeof(long)* np;
  memset(t->pevts.nums, 0, size);
  memset(t->pevts.iters, 0, size);
  memset(t->pevts.gids, 0, sizeof(long *)* np);
  memset(t->psends.nums, 0, size);
  memset(t->psends.iters, 0, size);
  memset(t->psends.gids, 0, sizeof(long *)* np);
  memset(t->precvs.nums, 0, size);
  memset(t->precvs.iters, 0, size);
  memset(t->precvs.gids, 0, sizeof(long *)* np);
  memset(t->pcolls.nums, 0, size);
  memset(t->pcolls.iters, 0, size);
  memset(t->pcolls.at, 0, sizeof(double (*)[2])* np);
  memset(t->pcolls.comm, 0, sizeof(long *)* np);
}

TraceData *CreateTrace(const ParaverFile *const file)
{
  TraceData *t= (TraceData *) malloc(sizeof(TraceData));
  memset(t, 0, sizeof(TraceData));

  t->runtime= ParaverFileGetRuntime(file);
  memcpy(t->timeunit, ParaverFileGetTimeUnit(file), 3);
  t->numnodes= ParaverFileGetNumNodes(file);
  t->numapps= ParaverFileGetNumApps(file);

  allocComms(t, ParaverFileGetNumComms(file), ParaverFileGetAllCommsSizes(file));
  ParaverFileReadComms(file, t->comms.sizes, t->comms.ranks);

  t->numprocs= ParaverFileGetNumProcs(file);

  allocLevel0Data(t);
  initLevel0Data(t);

  return t;
}

void allocLevel1Data()
{
  const int np= TraceGetNumProcs();
  ErrorIf(np< 1, "%s: Invlalid num-processes (%d) retrieved from trace header\n",
          __func__, np);

  /* evts */
  long num= 0;
  for(int ip= 0; ip< np; ++ip) {
    num+= TraceGetNumProcEvts(ip);
  }
  ErrorIf(num< 1, "%s: Invalid num-events (%ld) retrieved from trace body\n",
          __func__, num);
  TraceSetNumEvts(num);

  TraceSetPtrEvtsAt(malloc(sizeof(double)* num));
  TraceSetPtrEvtsId(malloc(sizeof(int)* num));

  TraceSetPtrEvtsCrit(malloc(sizeof(double)* num));
  TraceSetPtrEvtsProc(malloc(sizeof(int)* num));

  TraceSetPtrProcEvtsGids(Alloc2d_long(np, TraceGetPtrNumProcEvts(), num));

  /* sends */
  num= 0;
  for(int ip= 0; ip< np; ++ip) {
    num+= TraceGetNumProcSends(ip);
  }
  ErrorIf(num< 1, "%s: Invalid num-messages (%ld) retrieved from trace body\n",
          __func__, num);
  TraceSetNumMsgs(num);

  TraceSetPtrMsgsSendAt(malloc(sizeof(double[3])* num));
  TraceSetPtrMsgsRecvAt(malloc(sizeof(double[3])* num));
  TraceSetPtrMsgsSendRank(malloc(sizeof(int)* num));
  TraceSetPtrMsgsRecvRank(malloc(sizeof(int)* num));
  TraceSetPtrMsgsSize(malloc(sizeof(double)* num));
  TraceSetPtrMsgsTag(malloc(sizeof(int)* num));

  TraceSetPtrProcSendsGids(Alloc2d_long(np, TraceGetPtrNumProcSends(), num));

  /* recvs */
  TraceSetPtrProcRecvsGids(Alloc2d_long(np, TraceGetPtrNumProcRecvs(), num));

  /* pcolls */
  num= 0;
  for(int ip= 0; ip< np; ++ip) {
    num+= TraceGetNumProcColls(ip);
  }
  TraceSetNumAllProcColls(num);

  TraceSetPtrProcCollsAt(Alloc2d_double2(np, TraceGetPtrNumProcColls(), num));
  TraceSetPtrProcCollsComm(Alloc2d_int(np, TraceGetPtrNumProcColls(), num));
}

void initLevel1Data()
{
  long num= TraceGetNumEvts();
  memset(TraceGetPtrEvtsAt(), 0, sizeof(double)* num);
  memset(TraceGetPtrEvtsId(), 0, sizeof(int)* num);

  memset(TraceGetPtrEvtsCrit(), 0, sizeof(double)* num);
  /* for(long i= 0; i< num; ++i) { */
  /*   TraceSetEvtCrit(i, -1.0); */
  /* } */
  memset(TraceGetPtrEvtsProc(), 0, sizeof(int)* num);

  memset(TraceGetPtrProcEvtsGids(0), 0, sizeof(long)* num);

  num= TraceGetNumMsgs();
  memset(TraceGetPtrMsgsSendAt(), 0, sizeof(double[3])* num);
  memset(TraceGetPtrMsgsRecvAt(), 0, sizeof(double[3])* num);
  memset(TraceGetPtrMsgsSendRank(), 0, sizeof(int)* num);
  memset(TraceGetPtrMsgsRecvRank(), 0, sizeof(int)* num);
  memset(TraceGetPtrMsgsSize(), 0, sizeof(double)* num);
  memset(TraceGetPtrMsgsTag(), 0, sizeof(int)* num);

  memset(TraceGetPtrProcSendsGids(0), 0, sizeof(long)* num);
  memset(TraceGetPtrProcRecvsGids(0), 0, sizeof(long)* num);

  num= TraceGetNumAllProcColls();
  memset(TraceGetPtrProcCollsAt(0), 0, sizeof(double[2])* num);
  memset(TraceGetPtrProcCollsComm(0), 0, sizeof(int)* num);
}

void TraceAllocAndInitLevel1Data()
{
  allocLevel1Data();
  initLevel1Data();
}

TraceData *Trace0= NULL;

inline static void TraceConnectProcMsgIxToEvtMsgList(IndexList **const emList,
                                                     const long evtIter, IndexList *pmIx)
{
  IndexList *ptr= emList[evtIter];
  pmIx->next= ptr;
  emList[evtIter]= pmIx;
}
inline static void TraceConnectProcSendIxToEvtSends(const int p,
                                                    const double *const at, IndexList *ixs0, IndexList *ixs1)
{
  TraceConnectProcMsgIxToEvtMsgList(TraceGetPtrEvtSends(0),
                                    TraceGetGidProcEvtAt(p, at[0], -1), ixs0);
  TraceConnectProcMsgIxToEvtMsgList(TraceGetPtrEvtSends(1),
                                    TraceGetGidProcEvtAt(p, at[1], 1), ixs1);
}
inline static void TraceConnectProcRecvIxToEvtRecvs(const int p,
                                                    const double *const at, IndexList *ixr0, IndexList *ixr1)
{
  TraceConnectProcMsgIxToEvtMsgList(TraceGetPtrEvtRecvs(0),
                                    TraceGetGidProcEvtAt(p, at[0], -1), ixr0);
  TraceConnectProcMsgIxToEvtMsgList(TraceGetPtrEvtRecvs(1),
                                    TraceGetGidProcEvtAt(p, at[1], 1), ixr1);
}
void TraceConnectEvtsToMsgs()
{
  /* move to alloc-init-level-1 data later */

  const long npoints= TraceGetNumMsgs()* 2;
  Trace0->msgs.ipsmem= (IndexList *) malloc(sizeof(IndexList)* npoints);
  memset(Trace0->msgs.ipsmem, 0, sizeof(IndexList)* npoints);
  Trace0->msgs.iprmem= (IndexList *) malloc(sizeof(IndexList)* npoints);
  memset(Trace0->msgs.iprmem, 0, sizeof(IndexList)* npoints);

  const long nevts= TraceGetNumEvts();
  Trace0->evts.slist[0]= (IndexList **) malloc(sizeof(IndexList *)* nevts);
  memset(Trace0->evts.slist[0], 0, sizeof(IndexList *)* nevts);
  Trace0->evts.slist[1]= (IndexList **) malloc(sizeof(IndexList *)* nevts);
  memset(Trace0->evts.slist[1], 0, sizeof(IndexList *)* nevts);

  Trace0->evts.rlist[0]= (IndexList **) malloc(sizeof(IndexList *)* nevts);
  memset(Trace0->evts.rlist[0], 0, sizeof(IndexList *)* nevts);
  Trace0->evts.rlist[1]= (IndexList **) malloc(sizeof(IndexList *)* nevts);
  memset(Trace0->evts.rlist[1], 0, sizeof(IndexList *)* nevts);

  const int np= TraceGetNumProcs();

  if(GlOpts.show_opts.diag> 0) {
    long total= 0;
    for(int ip= 0; ip< np; ++ip) {
      total+= TraceGetNumProcSends(ip);
    }
    if(TraceGetNumMsgs()!= total) {
      Error("sum(#proc-sends) (%ld) mismatches with #msgs (%ld)\n", total,
            TraceGetNumMsgs());
    } else {
      Debug1("sum(#proc-sends) (%ld) matches with #msgs (%ld)\n", total,
             TraceGetNumMsgs());
    }

    total= 0;
    for(int ip= 0; ip< np; ++ip) {
      total+= TraceGetNumProcRecvs(ip);
    }
    if(TraceGetNumMsgs()!= total) {
      Error("sum(#proc-recvs) (%ld) mismatches with #msgs (%ld)\n", total,
            TraceGetNumMsgs());
    } else {
      Debug1("sum(#proc-recvs) (%ld) matches with #msgs (%ld)\n", total,
             TraceGetNumMsgs());
    }
  }

  long ixIt= 0;
  IndexList *ixl= TraceGetPtrMsgsProcSendIndex();
  for(int ip= 0; ip< np; ++ip) {
    for(long ips= 0; ips< TraceGetNumProcSends(ip); ++ips) {
      const double *const t= TraceGetProcSendAts(ip, ips);
      ixl[ixIt].i= ips; ixl[ixIt+ 1].i= ips;
      TraceConnectProcSendIxToEvtSends(ip, t, ixl+ ixIt, ixl+ ixIt+ 1);
      ixIt+= 2;
    }
  }
  ErrorIf(TraceGetNumMsgs()* 2!= ixIt,
          "#messages= %ld, iterator(x2)= %ld\n", TraceGetNumMsgs(), ixIt);

  ixIt= 0;
  ixl= TraceGetPtrMsgsProcRecvIndex();
  for(int ip= 0; ip< np; ++ip) {
    for(long ipr= 0; ipr< TraceGetNumProcRecvs(ip); ++ipr) {
      double *t= TraceGetProcRecvAts(ip, ipr);
      ixl[ixIt].i= ipr; ixl[ixIt+ 1].i= ipr;
      TraceConnectProcRecvIxToEvtRecvs(ip, t, ixl+ ixIt, ixl+ ixIt+ 1);
      ixIt+= 2;
    }
  }
  ErrorIf(TraceGetNumMsgs()* 2!= ixIt,
          "#messages= %ld, iterator(x2)= %ld\n", TraceGetNumMsgs(), ixIt);
}
