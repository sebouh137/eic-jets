reset
set terminal epslatex size 20cm,6.75cm color dashed dl 3 standalone
LEFT= 0.15
DX=0.25
set style fill  transparent solid 0.50 noborder
set output 'positron_beam_with_stats.tex'
set multiplot layout 1,3

#------1st plot-----
set lmargin at screen LEFT
set rmargin at screen LEFT+DX

set xtics nomirror
set ytics nomirror
set xtics 0.1
set mxtics 1
set mytics 1
set xrange [0.01:0.59]
set yrange [-0.125:0.125]
set style fill pattern border lc "black"
#set xlabel '$z$' offset 1,1
set xlabel '$z_h$' offset 0,0.2
set ylabel '$A_{UT}^{\sin(\phi_{S_A}-\hat{\phi}_h)}$' offset 1,0
#set ylabel '$xf(x,Q^2=10$GeV$^2)$' offset 1,0
set style line 1 lc rgb "red"   lw 4.0 lt 2 # red solid
set style line 2 lc rgb "black" lw 2.0 lt 0 # black dashed
set style line 6 lc rgb "black" lw 5.5 lt 0 # black dashed
set style line 3 lc rgb "black"  lw 4.0 dashtype 3 lt 1 # blue solid
set style line 4 lc rgb "violet" lw 2.0 dashtype 2 lt 2 # green dashed
set style line 5 lc rgb "blue" lw 4.0 dashtype 2 lt 4 # black dotted
set label '\scriptsize $\mathrm{\sqrt{s} = 13~TeV,~anti\mbox{-}k_T,~R=0.8}$' at -6.6,0.7875
set label '\scriptsize $\mathrm{p_T > 600~GeV,~|\eta|<1.5}$' at -6.0,0.7125
set label 'Positron beam' at 0.05,0.07
set label '$0.10<x<0.20$' at 0.05,0.045
#set label '\scriptsize CC DIS 10+275 GeV, 100 fb$^{-1}$, $0.01<y<0.9,\ Q^2>100$ GeV$^2$' at 0.35,0.16

#set key title '\scriptsize u quark, $Q^2=10$GeV$^2$'
#set key title '\scriptsize $0.10<x<0.20$' at 0.4,0.17

#set table 'lower1.dat'
#plot '../fig3_pip_e+_p1.dat' using 1:($2-$3) smooth sbezier
#unset table
#set table 'upper1.dat'
#plot '../fig3_pip_e+_p1.dat' using 1:($2+$3) smooth sbezier
#unset table

#set table 'lower2.dat'
#plot '../fig3_pim_e+_p1.dat' using 1:($2-$3) smooth sbezier
#unset table
#set table 'upper2.dat'
#plot '../fig3_pim_e+_p1.dat' using 1:($2+$3) smooth sbezier
#unset table

#set table 'lower3.dat'
#plot '../fig3_pip_e+_p2.dat' using 1:($2-$3) smooth sbezier
#unset table
#set table 'upper3.dat'
#plot '../fig3_pip_e+_p2.dat' using 1:($2+$3) smooth sbezier
#unset table

#set table 'lower4.dat'
#plot '../fig3_pim_e+_p2.dat' using 1:($2-$3) smooth sbezier
#unset table
#set table 'upper4.dat'
#plot '../fig3_pim_e+_p2.dat' using 1:($2+$3) smooth sbezier
#unset table

#set table 'lower5.dat'
#plot '../fig3_pip_e+_p3.dat' using 1:($2-$3) smooth sbezier
#unset table
#set table 'upper5.dat'
#plot '../fig3_pip_e+_p3.dat' using 1:($2+$3) smooth sbezier
#unset table

#set table 'lower6.dat'
#plot '../fig3_pim_e+_p3.dat' using 1:($2-$3) smooth sbezier
#unset table
#set table 'upper6.dat'
#plot '../fig3_pim_e+_p3.dat' using 1:($2+$3) smooth sbezier
#unset table

plot '../stats/211_e+_p0.dat' u 1:(0):3 with errorbars t '' pt 6 ps 1.5 lw 2.5 lt rgb "black" ,\
     '../stats/-211_e+_p0.dat' u ($1+.01):(0):3 with errorbars t '' pt 7 ps 1.5 lw 2.5 lt rgb "red" ,\
     'y0.dat' u 1:2 w lines lc rgb "blue"  lw 5 title '' ,\

unset label
#------2nd plot-----
set ytics format ""
#set key title '\scriptsize d quark, $Q^2=10$GeV$^2$'
unset ylabel
set lmargin at screen LEFT+1*DX
set rmargin at screen LEFT+2*DX

set label '\scriptsize Single inclusive groomed jet' at -6.4,0.875
set label '\scriptsize $\mathrm{\beta = 1},~a=-1$' at -5.0,0.625
#set key title '\scriptsize $0.20<x<0.30$'
set label '$0.20<x<0.30$' at 0.05,0.045

plot '../stats/211_e+_p1.dat' u 1:(0):3 with errorbars t '$\pi^+$' pt 6 ps 1.5 lw 2.5 lt rgb "black" ,\
     '../stats/-211_e+_p1.dat' u ($1+.01):(0):3 with errorbars t '$\pi^-$' pt 7 ps 1.5 lw 2.5 lt rgb "red" ,\
     'y0.dat' u 1:2 w lines lc rgb "blue" lw 5 title 'prediction' ,\


unset ylabel
unset label

#------3nd plot-----
set ytics format ""
#set key title '\scriptsize $0.30<x<0.50$'
unset ylabel
set lmargin at screen LEFT+2*DX
set rmargin at screen LEFT+3*DX

set label '\scriptsize Single inclusive groomed jet' at -6.4,0.875
set label '\scriptsize $\mathrm{\beta = 1},~a=-1$' at -5.0,0.625

set label '$0.30<x<0.50$' at 0.05,0.045

plot '../stats/211_e+_p2.dat' u 1:(0):3 with errorbars t '' pt 6 ps 1.5 lw 2.5 lt rgb "black" ,\
     '../stats/-211_e+_p2.dat'  u ($1+.01):(0):3 with errorbars t '' pt 7 ps 1.5 lw 2.5 lt rgb "red" ,\
     'y0.dat' u 1:2 w lines lc rgb "blue" lw 5 title '' ,\


unset ylabel
unset label


unset multiplot
set output # finish the current output file
system('latex positron_beam_with_stats.tex && dvips positron_beam_with_stats.dvi && ps2eps -f positron_beam_with_stats.ps && epstopdf positron_beam_with_stats.eps && rm positron_beam_with_stats-inc.eps positron_beam_with_stats.ps *.aux *.dvi *.log *.tex')
