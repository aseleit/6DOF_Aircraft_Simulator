set datafile separator ","
set autoscale fix
set key outside right center
set grid

# First plot
set terminal qt 0
set xlabel 'Time (s)' font "Times-Roman,14"
set ylabel 'Pitch Rate (deg/s)' font "Times-Roman,14"

plot 'sim_out.csv' using 1:2 with lines lc rgb "red" lw 2

pause -1