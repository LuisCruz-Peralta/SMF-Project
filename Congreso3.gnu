set terminal pngcairo size 768, 576 font "Liberation Serif Bold, 20"

stats 'Congreso.txt' using 1:2
x_lower = STATS_min_x
x_upper = STATS_max_x

y_lower = STATS_min_y
y_upper = STATS_max_y

stats 'Congreso.txt' using 3
z_lower = STATS_min
z_upper = STATS_max

stop_run = STATS_records

i = 1

do for [k = 1:stop_run:5] {
set output sprintf('./Outputs/lorenz_gnuplot_%d.png',i)
i = i + 1
set title "Simulacion"
set tics font "Liberation Serif, 12"
set xrange [x_lower:x_upper]
set yrange [y_lower:y_upper]
set zrange [z_lower:z_upper]
set view 75,10
splot 'Congreso.txt' every ::1::k with lines linecolor "midnight-blue" title ""
unset output
}