set datafile separator ","
set autoscale fix
set key outside right center
set grid
set encoding utf8

# First plot
set terminal qt 0
set xlabel 'Time (s)' font "Times-Roman,14"

plot 'sim_out.csv' using 1:2 with lines lc rgb "red" lw 2 title "{/Symbol D} u (m/s)",\
     'sim_out.csv' using 1:6 with lines lc rgb "green" lw 2 title "{/Symbol D} v (m/s)",\
     'sim_out.csv' using 1:3 with lines lc rgb "black" lw 2 title "{/Symbol D} w (m/s)"


# Second plot
set terminal qt 1
set xlabel 'Time (s)' font "Times-Roman,14"

plot 'sim_out.csv' using 1:7 with lines lc rgb "red" lw 2 title "{/Symbol D} p (deg/s)",\
     'sim_out.csv' using 1:4 with lines lc rgb "green" lw 2 title "{/Symbol D} q (deg/s)",\
     'sim_out.csv' using 1:8 with lines lc rgb "black" lw 2 title "{/Symbol D} r (deg/s)"


# Third plot
set terminal qt 2
set xlabel 'Time (s)' font "Times-Roman,14"

plot 'sim_out.csv' using 1:9 with lines lc rgb "red" lw 2 title "{/Symbol D} {/Symbol f} (deg)",\
     'sim_out.csv' using 1:5 with lines lc rgb "green" lw 2 title "{/Symbol D} {/Symbol q} (deg)",\
     'sim_out.csv' using 1:10 with lines lc rgb "black" lw 2 title "{/Symbol D} {/Symbol y} (deg)"

pause -1