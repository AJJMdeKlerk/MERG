set terminal epslatex color standalone dash
set out 'plot.tex'
set key Left
set key width 3
set key height 1
set key spacing 1.25
set key top right
set key reverse
set key font 'Times-Roman,18'
set xrange [0:50000]
set yrange [0.0:0.5]
set ytics 
set xlabel 'Number of basis states' font 'Times-Roman, 30'
set ylabel '$(E_0-E_{\mathrm{exact}})/E_{\mathrm{exact}}$'
#set y2label 'Two particle weight, $w_2$' offset -2
set xtics  font 'Times-Roman, 22'
set ytics  font 'Times-Roman, 22'
#set y2tics font 'Times-Roman, 22'
set tics scale 1.5
set mytics
set format x "%2.1t\$\\times10^{%T}$"
#set size 1, 2
set for [i=1:5] xtics (0,10000*i)
set xtics add ('0' 0)
#set mxtics

#set datafile separator ", "
# plot 'infile' using 0:1

plot "./Energies_statesKept_720_statesAdded_80_maxNRGSteps_615_method_HOSTS.csv" u ($1):($2) w lp lw 1 lc 4 pt 4 t 'HOSTS', "./Full.csv" u ($1):($2) w p pt 7 lc -1 t 'Full' 
