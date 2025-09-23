#
# Copyright (c) 2025      High Performance Computing Center Stuttgart,
#                         University of Stuttgart.  All rights reserved.
#
# Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
#

set title '{/:Bold Local and global metrics}' font ',18'
set xlabel '{/:Bold Elapsed time [s]}' font ',18'
set ylabel '{/:Bold Metrics value}' font ',18' offset -1,0
set grid
set key default
set key bottom right Right width -7 font ',24'

set xtics font ',18'
set ytics font ',18'
set lmargin 10
set bmargin 4

p fname u ($1*1.0e-9):($5/$4) w lp lw 4 lc rgb 'dark-cyan' title 'LB^@{}_{local}'\
, '' u ($1*1.0e-9):($4/$7) w lp lw 4 lc rgb 'forest-green' title 'Ser^@{}_{local}'\
, '' u ($1*1.0e-9):($7/$6) w lp lw 4 lc rgb 'coral' title 'Trf^@{}_{local}'

#, <global-lb-val> w l lc rgb 'dark-cyan' dt 2 lw 2 title 'LB_{global}', <global-ser-val> w l lc rgb 'forest-green' dt 2 lw 2 title 'Ser_{global}', <global-trf-val> w l lc rgb 'coral' dt 2 lw 2 title 'Trf_{global}'

ROUND_X=.03*(GPVAL_DATA_X_MAX- GPVAL_DATA_X_MIN)
set xrange[GPVAL_DATA_X_MIN-ROUND_X:GPVAL_DATA_X_MAX+ROUND_X]
set yrange[-0.05:1.05]

#uncomment for png output; do not need to use persist in that case
#set term pngcairo enhanced dashed size 800,600
#set output fname[:strlen(fname)-4].'.png'

replot
