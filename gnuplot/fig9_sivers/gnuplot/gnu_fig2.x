reset
set terminal epslatex size 50cm,10cm color dashed dl 3 standalone font "cmr,11"
LEFT= 0.15
DX=0.175
set style fill  transparent solid 0.50 noborder
set output 'sivers_theory_with_stats.tex'
set multiplot layout 1,1

#------1st plot-----
set lmargin at screen LEFT
set rmargin at screen LEFT+DX

set xtics nomirror
set ytics nomirror
set xtics 
set ytics 
set mxtics 1
set mytics 1
set xrange [0.05:0.85]
set yrange [-6.2:8]
set style fill pattern border lc "black"
set xlabel '$x$' offset 0,0.2
set ylabel '$A_{UT}^{\sin(\phi_q-\phi_{S_A})}\ [\%]$' offset 0,0
set style line 1 lc rgb "red"   lw 3.0 lt 2 # red solid
set style line 2 lc rgb "black" lw 2.0 lt 0 # black dashed
set style line 6 lc rgb "black" lw 5.5 lt 0 # black dashed
set style line 3 lc rgb "black"  lw 3.0 lt 1 # blue solid
set style line 4 lc rgb "dark-violet" lw 3.0 lt 4 # black dotted
set style line 5 lc rgb "blue" lw 3.0 lt 4 # black dotted
#set label '\scriptsize $Q^2>100$ GeV$^2$' at 0.45,-1.7
#set label '\scriptsize $0.1<y<0.9$' at 0.45,-2.5
#set label '\scriptsize $0<q_T<5$ GeV' at 0.45,-3.3
#set label '\scriptsize $\sqrt{s}=105$ GeV' at 0.45,-4.1
#set label '\scriptsize $R=1.0$' at 0.45,-4.9
set label '$Q^2>100$ GeV$^2$' at 0.25,-3.9
set label '$0.1<y<0.9$' at 0.25,-4.8
set label '$0<q_T<5$ GeV' at 0.25,-5.7
set label '$\sqrt{s}=105$ GeV' at 0.57,-4.8
set label '$R=1.0$' at 0.57,-5.7


set table 'lower1.dat'
plot '../fig2_w+.dat' using 1:(($2+$3)*100) smooth csplines
unset table
set table 'upper1.dat'
plot '../fig2_w+.dat' using 1:(($2-$3)*100) smooth csplines
unset table

set table 'lower2.dat'
plot '../fig2_w-.dat' using 1:(($2+$3)*100) smooth csplines
unset table
set table 'upper2.dat'
plot '../fig2_w-.dat' using 1:(($2-$3)*100) smooth csplines
unset table

plot '< paste lower1.dat upper1.dat' using 1:2:5 with filledcurves fs solid 0.4 lt rgb "orange" title '' ,\
     '< paste lower2.dat upper2.dat' using 1:2:5 with filledcurves fs solid 0.4 lt rgb "blue" title '' ,\
      '../fig2_w+.dat' u 1:($2*100) w lines smooth csplines lc 'orange' lw 5 title 'Positron beam',\
      'stats_positron.dat' with errorbars title 'Stat. precision' pt 6 ps 1.5 lw 2.5 lt rgb "black",\
      '../fig2_w-.dat' u 1:($2*100) w lines smooth csplines lc 'blue' lw 5 title 'Electron beam',\
      'stats_electron.dat' using ($1+0.015):2:3 with errorbars title 'Stat. precision' pt 7 ps 1.5 lw 2.5 lt rgb "black",\
      'y0.dat' w lines lw 4 lt '-' lc 'black' title ''

#plot 'stats_electron.dat' with errorbars title '' pt 7 ps 2 lw 3 lt rgb "black"
#plot 'stats_positron.dat' using ($1+0.015):2:3 with errorbars title '' pt 6 ps 2 lw 3 lt rgb "black"

unset ylabel
unset label
set ytics format ""

unset multiplot
set output # finish the current output file
system('latex sivers_theory_with_stats.tex && dvips sivers_theory_with_stats.dvi && ps2eps -f sivers_theory_with_stats.ps && epstopdf sivers_theory_with_stats.eps && rm sivers_theory_with_stats.inc.eps sivers_theory_with_stats.ps *.aux *.dvi *.log *.tex')
