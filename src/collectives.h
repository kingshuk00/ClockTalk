/*
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef REPLAY_COLLECTIVES_H__
#define REPLAY_COLLECTIVES_H__

#include"common.h"
#include"trace_data.h"
#include"paraver.h"
#include"clocks.h"
#include<float.h>

static struct {
  double *entry;   /* critical time of entry */
  bool *pexec;     /* whether proc is part of the comm - constant*/
  double last;
  int nremains;    /* #memebers remains to enter in this collective */
  int evt;         /* collective-event-id */
} *Colls= NULL;

static void CollsAlloc()
{
  const long nc= TraceGetNumComms();
  Colls= malloc(sizeof(*Colls)* nc);
  memset(Colls, 0, sizeof(*Colls)* nc);

  long ncvalid= 0;
  for(long c= 0; c< nc; ++c) {
    if(TraceIsCommSelf(c)) {
      continue;  /* no need for counting COMM_SELFs */
    }
    ++ncvalid;
  }
  const int np= TraceGetNumProcs();
  double *entry= (double *) malloc(sizeof(double)* ncvalid* np);
  memset(entry, 0, sizeof(double)* ncvalid* np);
  bool *pexec= (bool *) malloc(sizeof(bool)* ncvalid* np);
  memset(pexec, 0, sizeof(bool)* ncvalid* np);
  double *eptr= entry;
  bool *xptr= pexec;
  for(long c= 0; c< nc; ++c) {
    Colls[c].evt= -1;
    Colls[c].last= -1.0;

    if(TraceIsCommSelf(c)) {
      continue;  /* COMM_SELF */
    }

    Colls[c].entry= eptr;
    eptr+= np;

    for(int i= 0; i< TraceGetCommSize(c); ++i) {
      xptr[TraceGetCommRank(c, i)]= true;
    }
    Colls[c].pexec= xptr;
    xptr+= np;
  }
}
#if 0
static void CollsFree()
{
  if(NULL== Colls) {
    return;
  }
  FREE_IF(Colls[0].executing);
  FREE_IF(Colls[0].critical);
  FREE_IF(Colls[0].elapsed);
  FREE_IF(Colls);
}
#endif

static void CollsResetOneColl(const int c)
{
  Debug1("Resetting collective of comm-%d\n", c);
  Colls[c].nremains= TraceGetCommSize(c);
  Colls[c].evt= -1;
  Colls[c].last= -1.0;
  if(TraceIsCommSelf(c)) {
    return;
  }
  for(int ip= 0; ip< TraceGetCommSize(c); ++ip) {
    Colls[c].entry[TraceGetCommRank(c, ip)]= -1.0;
  }
}
static void CollsResetAll()
{
  for(int c= 0; c< TraceGetNumComms(); ++c) {
    CollsResetOneColl(c);
  }
}

inline static int GetCollEvt(const int c) { return Colls[c].evt; }
inline static bool IsCollActive(const int c) { return Colls[c].evt> 0; }
inline static bool IsCollAvailable(const int c, const int p) { return -1.0!= Colls[c].entry[p]; }
inline static const char *GetCollName(const int c) { return GetParaverMPIEvtName(Colls[c].evt); }
inline static void ActivateColl(const int c, const int p, const int evt)
{
  ErrorIf(TraceIsCommSelf(c), "%d: activating coll %s with COMM_SELF!\n", p,
          GetParaverMPIEvtName(evt));
  ErrorIf(IsCollActive(c),
          "%d: activating coll %s: already active collective %s(%d)\n",
          p, GetParaverMPIEvtName(evt), GetCollName(c), c);

  Debug1("%d: trigger coll %s(%d) - %d/%d\n", p, GetParaverMPIEvtName(evt), c,
         Colls[c].nremains, TraceGetCommSize(c));
  Colls[c].evt= evt;
}

inline static bool EveryoneEnteredColl(const int c) { return 0== Colls[c].nremains&& TraceGetCommSize(c)> 0; }
inline static bool LastCollEntryEstablished(const int c) { return -1.0!= Colls[c].last; }

inline static void EnterColl(const int c, const int p, const double t,
                             const int evt)
{
  ErrorIf(-1.0!= Colls[c].entry[p],
          "%d: coll %s(%d) since crit-%.0lf, overwrite at crit-%.0lf (elapsed-%.0lf)\n",
          p, GetCollName(c), c, Colls[c].entry[p], ClockGetCritical(p), t);

  if(!IsCollActive(c)) {
    ActivateColl(c, p, evt);
  }
  Colls[c].entry[p]= ClockGetCritical(p);
  --(Colls[c].nremains);

  Debug1("%d: coll %s(%d) enter at %.0lf (critical: %.0lf) - %d/%d remains\n", p,
         GetCollName(c), c, t, ClockGetCritical(p), Colls[c].nremains,
         TraceGetCommSize(c));

  if(EveryoneEnteredColl(c)) {
    ErrorIf(-1.0!= Colls[c].last,
            "%d: everyone just entered coll %s(%d), but last entry already set at %.0lf\n",
            p, GetCollName(c), c, Colls[c].last);
    double last= 0.0;
    for(int ip= 0; ip< TraceGetCommSize(c); ++ip) {
      last= MAX(Colls[c].entry[TraceGetCommRank(c, ip)], last);
    }
    Debug1("%d: last-critical-entry into coll %s at %.0lf\n", p, GetCollName(c),
           last);
    Colls[c].last= last;
  }
}
/* return 0 if finished, 1 otherwise */
inline static int LeaveColl(const int c, const int p, const double t,
                            const int collEvt)
{
  if(!LastCollEntryEstablished(c)) {
    return 1;
  }
  ClockUpdateCritical(p, Colls[c].last- Colls[c].entry[p]);
  ++(Colls[c].nremains);
  Debug1("%d: coll %s(%d) leave at %.0lf (critical: %.0lf) - %d/%d done\n", p,
         GetCollName(c), c, t, ClockGetCritical(p), Colls[c].nremains,
         TraceGetCommSize(c));
  if(TraceGetCommSize(c)== Colls[c].nremains) {
    CollsResetOneColl(c);
  }
  return 0;
}

#endif  /* REPLAY_COLLECTIVES_H__ */
