reset


set style line 1 lc rgb "#ff0000" lt 1 lw 1.5
set style line 2 lc rgb "#ffff00" lt 1 lw 1.5
set style line 3 lc rgb "#00ff00" lt 1 lw 1.5
set style line 4 lc rgb "#0000FF" 
set style line 5 lc rgb "#ff0000" lw 3
set style line 6 lc rgb "#ff0000" lt 2 lw 3



  
set xrange [0:8]
set yrange [0:1]
set xlabel "r"
set ylabel "phi"

	plot "phi.dat" using 1:2 ls 1 title "phiA1", \
"phi.dat" using 1:3 ls 2 title "phiB1", \
"phi.dat" using 1:4 ls 5 title "phiA2", \
"phi.dat" using 1:5 ls 3 title "phiB2", \
"phi.dat" using 1:6 ls 6 title "phiA3", \
"phi.dat" using 1:7 ls 4 title "phiC", \
