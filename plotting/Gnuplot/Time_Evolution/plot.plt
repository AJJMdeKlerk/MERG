set terminal epslatex color standalone dash
set out 'plot.tex'
set key Left
set key width 3
set key height 1
set key spacing 1.25
set key bottom right
set key reverse
set key font 'Times-Roman,18'
set logscale xy
set xrange [0.001:1]
set yrange [0.001:1]
set ytics 
set xlabel 'Time $k_F^2t$' font 'Times-Roman, 30'
set ylabel '$g_2(t)$'
#set y2label 'Two particle weight, $w_2$' offset -2
set xtics  font 'Times-Roman, 22'
set ytics  font 'Times-Roman, 22'
#set y2tics font 'Times-Roman, 22'
set tics scale 1.5
set mytics
#set format x "%2.1t\$\\times10^{%T}$"
#set size 1, 2
#set for [i=1:5] xtics (0,0.2*i)
set xtics add ('0' 0)
#set mxtics

#set datafile separator ", "
# plot 'infile' using 0:1

plot "./NJP_ZillEtAl_Extracted_N5_data.txt" u ($1 / 6.3165468167):($2) w l lw 6 lc 7 t 'Bethe Ansatz', "./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_10000.csv" u ($1):($2) w l lw 6 lc 4 t 'MERG $N_{tot}$ = 10000', "./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_20000.csv" u ($1):($2) w l lw 6 lc 5 t 'MERG $N_{tot}$ = 20000', "./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_30000.csv" u ($1):($2) w l lw 6 lc 2 t 'MERG $N_{tot}$ = 30000', "./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_40000.csv" u ($1):($2) w l lw 6 lc 6 t 'MERG $N_{tot}$ = 40000', "./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_50000.csv" u ($1):($2) w l lw 6 lc 9 t 'MERG $N_{tot}$ = 50000' 

#"./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_5000.csv" u ($1):($2) w l lw 6 lc 5 t 'MERG $N_{tot}$ = 5000',
#"./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_15000.csv" u ($1):($2) w l lw 6 lc 7 t 'MERG $N_{tot}$ = 15000',
#"./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_25000.csv" u ($1):($2) w l lw 6 lc 9 t 'MERG $N_{tot}$ = 25000',
#"./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_35000.csv" u ($1):($2) w l lw 6 lc 1 t 'MERG $N_{tot}$ = 35000',
#"./TimeEvolution_totalTime_1_timeStepSize_0.0001__statesKept700_statesAdded_100_maxNRGSteps_492_method_MERG_currentBasisSize_45000.csv" u ($1):($2) w l lw 6 lc -1 t 'MERG $N_{tot}$ = 45000',
