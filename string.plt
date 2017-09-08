#!/usr/bin/gnuplot

set terminal postscript color enhanced
set size 0.7, 0.7
set xrange [-1:11]
set yrange [-1:11]
set zrange [-5:6]

set output "string.eps"
set ticslevel 0
splot "string.dat" u 1:2:3 with lines