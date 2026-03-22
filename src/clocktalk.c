/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */


#include"common.h"
#include"build_info.h"
#include"utils/utils.h"
#include"opt_parser.h"
#include"paraver/paraver.h"
#include"replay/replay.h"
#include"trace/trace.h"
#include"clocks/clocks.h"
#include"monitoring/monitoring.h"

int main(int argc, char *argv[])
{
  ClockTalkOpts *opts= ParseOpts(argc, argv);

  Debug1("Running program built on %s at %s\n", CT_BUILD_DATE, CT_BUILD_TIME);

  const double t0= Timer_s();
  if(0!= ReadParaverFile(opts)) {
    Error("Problem reading paraver file \"%s\"\n", argv[1]);
    return 0;
  }
  const double t1= Timer_s();

  if(opts->show.timings) {
    printf("Reading Paraver file took %.1lf s\n", t1- t0);
  }

  ReplayTrace(opts);
  if(false) {
    FILE *fp= fopen("checking.txt", "w");
    for(TraceResetIterEvts(); TraceGetIterEvts()< TraceGetNumEvts();
        TraceIncrIterEvts()) {
      fprintf(fp, "%.9e %.9e %3d\n",
              TraceGetCurrEvtAt(), TraceGetCurrEvtCrit(), TraceGetCurrEvtId());
    }
    fclose(fp); fp= NULL;
  }

  ShowStatsPostReplay(opts->show.pretty);

  ClockFinalize();

  if(opts->show.timings) {
    const double t2= Timer_s();
    printf("Replay took %.1lf s (total %.1lf s)\n", t2- t1, t2- t0);
  }

  if(opts->mon.evt.isOn) {
    DoMonitoringEventBased(opts);
  }

  if(opts->mon.win.isOn) {
    DoMonitoringWindowed(opts);
  }

  FREE_IF(opts->filename);
  FREE_IF(opts);

  return 0;
}
