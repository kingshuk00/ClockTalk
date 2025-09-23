/*
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef REPLAY_PARAVER_H__
#define REPLAY_PARAVER_H__

#include<stdbool.h>

extern int ReadParaverFile(const char *const);
extern const char *GetParaverMPIEvtName(const int);

inline static bool ParaverEvtIsColl(const int evtId)
{
  bool ret= false;
  switch(evtId) {
  case 7:   /* Bcast */
  case 8:   /* Barrier */
  case 9:   /* Reduce */
  case 10:  /* Allreduce */
  case 11:  /* Alltoall */
  case 12:  /* Alltoallv */
  case 13:  /* Gather */
  case 14:  /* Gatherv */
  case 17:  /* Allgather */
  case 18:  /* Allgatherv */
  case 30:  /* Scan */
  case 80:  /* Reduce_scatter */
  case 163: /* MPI_Igatherv */
    ret= true;
    break;
  default:
    break;
  }
  return ret;
}

/* Dimemas treats the turned off functions as 'local' - i.e. finishes instantly */
inline static bool ParaverCollEvtIsDimemasCompliant(const int evtId)
{
  bool ret= false;
  switch(evtId) {
  case 7:   /* Bcast */
  case 8:   /* Barrier */
  case 9:   /* Reduce */
  case 10:  /* Allreduce */
  case 11:  /* Alltoall */
  case 12:  /* Alltoallv */
  case 13:  /* Gather */
  case 14:  /* Gatherv */
  case 17:  /* Allgather */
#if 0  /* Dimemas compliance */
  case 18:  /* Allgatherv */
  case 30:  /* Scan */
#endif
  case 80:  /* Reduce_scatter */
  case 163: /* MPI_Igatherv */
    ret= true;
    break;
  default:
    break;
  }
  return ret;
}

#endif  /* REPLAY_PARAVER_H__ */
