#!/usr/bin/gnuplot -persist

set title "Simulating bodies"
set xlabel "Threads"
set ylabel "Execution time"
set grid
set terminal qt 0
plot "./benchmark/4_80000.dat" with lines, \
     "./benchmark/4_110000.dat" with lines, \
     "./benchmark/4_140000.dat" with lines, \
     "./benchmark/4_170000.dat" with lines, \
     "./benchmark/4_200000.dat" with lines, \
     "./benchmark/4_230000.dat" with lines
