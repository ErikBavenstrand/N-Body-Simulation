#!/usr/bin/gnuplot -persist

set title "Simulating bodies"
set xlabel "Threads"
set ylabel "Execution time"
set grid
set terminal qt 0
plot "./benchmark/2_250.dat" with lines, \
     "./benchmark/2_300.dat" with lines, \
     "./benchmark/2_350.dat" with lines, \
     "./benchmark/2_400.dat" with lines, \
     "./benchmark/2_450.dat" with lines, \
     "./benchmark/2_500.dat" with lines, \