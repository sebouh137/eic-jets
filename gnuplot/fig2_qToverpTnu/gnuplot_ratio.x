reset
set terminal epslatex size 50cm,10cm color dashed dl 3 standalone font "cmr,17"
LEFT= 0.15
DX=0.25
set style fill  transparent solid 0.50 noborder 
set output 'qToverptnu_theory.tex'
set multiplot layout 1,1

#------1st plot-----
set lmargin at screen LEFT
set rmargin at screen LEFT+DX

set xtics nomirror
set ytics nomirror
set xtics 
set ytics 
set mxtics 5
set mytics 5
set xrange [0.:0.4]
set yrange [0:9.5]
set style fill pattern border lc "black"
set xlabel '$q_T/p_T^{\nu}$' offset 0,0.2
set ylabel '$1/\sigma\times d\sigma/d[q_T/p_T^{\nu}]$' offset 0,0
set style line 1 lc rgb "red"   lw 3.0 lt 2 # red solid
set style line 2 lc rgb "black" lw 2.0 lt 0 # black dashed
set style line 6 lc rgb "black" lw 5.5 lt 0 # black dashed
set style line 3 lc rgb "black"  lw 3.0 lt 1 # blue solid
set style line 4 lc rgb "dark-violet" lw 3.0 lt 4 # black dotted
set style line 5 lc rgb "blue" lw 3.0 lt 4 # black dotted

set table 'up_smooth'
plot 'data/ratio_theory_new.dat' using 1:6 smooth csplines

set table 'down1_smooth'
plot 'data/ratio_theory_new.dat' using 1:4 smooth csplines

set table 'down2_smooth'
plot 'data/ratio_theory_new.dat' using 1:5 smooth csplines
unset table

plot '< paste up_smooth down1_smooth' u 1:2:5 with filledcurves fillstyle solid 0.5 lw 4.5 lc rgb 'orange' title '',\
     '< paste up_smooth down2_smooth' u 1:2:5 with filledcurves fillstyle solid 0.5 lw 4.5 lc rgb 'orange' title '',\
     'data/pythia_ratio.dat' w lines lc 'blue' lw 5 title '$\textsc{Pythia8}$',\
     'data/ratio_theory_new.dat' u 1:2 w lines smooth csplines lc 'orange' lw 5 title 'NLL',\
#     'fig5.dat' u 1:($2/$3) w lines smooth csplines lc 'orange' lw 5 title 'theory',\
#     'fig5_ave.dat' u 1:2 lc 'red' lw 5 title '',\
#     'y0.dat' w lines lc 'black' title ''
#     './j3_R1_frac.dat' u 1:($2*100) smooth csplines lt 7 lw 5 title '$0.5 < z_\Lambda < 0.8$',\
#     './0125_zh_5to8.dat' u 1:2 smooth csplines lt 6 lw 5 title '$0.5 < z_h < 0.8$',\
#     './zh_1to5.dat' u 1:2 w points ls 3 title '',\
#     './zh_5to8.dat' u 1:2 w points ls 4 title '',\
#     './teq_zh_1to5.dat' u 1:2 notitle,\
#     './cteq_zh_5to8.dat' u 1:2 notitle,\

unset ylabel
unset label
set ytics format ""

unset multiplot
set output # finish the current output file
system('latex qToverptnu_theory.tex && dvips qToverptnu_theory.dvi && ps2eps -f qToverptnu_theory.ps && epstopdf qToverptnu_theory.eps && rm qToverptnu_theory.inc.eps qToverptnu_theory.ps *.aux *.dvi *.log *.tex')
