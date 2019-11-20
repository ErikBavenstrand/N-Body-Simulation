#!/usr/bin/gnuplot -persist

set title "Simulating bodies"
set xlabel "Bodies"
set ylabel "Execution time"
set grid
set terminal qt 0
plot "./benchmark/1.dat" with lines