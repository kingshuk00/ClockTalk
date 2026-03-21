/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */

#include"utils.h"
#include<stdio.h>

static int upNo(const char *restrict format, ...) { return 0; }

int (*upErr)(const char *restrict format, ...)= upNo;
int (*upDbg1)(const char *restrict format, ...)= printf;
int (*upDbg2)(const char *restrict format, ...)= printf;
int (*upDbg3)(const char *restrict format, ...)= printf;
int (*upDbg4)(const char *restrict format, ...)= printf;
int (*upDbg5)(const char *restrict format, ...)= printf;

int UtilSetShowFunctions(const int errLevel, const int dbgLevel)
{
  if(errLevel> 0) {
    upErr= printf;
  }

  switch(dbgLevel) {
  case 0:
    upDbg1= upNo;               /* fall through */
  case 1:
    upDbg2= upNo;               /* fall through */
  case 2:
    upDbg3= upNo;               /* fall through */
  case 3:
    upDbg4= upNo;               /* fall through */
  default:
    upDbg5= upNo;               /* fall through */
    break;
  }

  return 0;
}
