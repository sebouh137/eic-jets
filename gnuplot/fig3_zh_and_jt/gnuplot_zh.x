reset
set terminal epslatex size 50cm,10cm color dashed dl 3 standalone font "cmr,14"
LEFT= 0.15
DX=0.175
set style fill  transparent solid 0.50 noborder
set output 'zh_theory.tex'
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
set logscale x
set logscale y
set xrange [0.05:1]
set yrange [0.001:35]
set style fill pattern border lc "black"
set xlabel '$z_h=(\vec{p}_{\textrm{jet}}\cdot\vec{p}_{\textrm{hadron}})/|\vec{p}_{\textrm{jet}}|^2$' offset 0,0.2
set ylabel '$1/\sigma\times d\sigma/dz_h$' offset 2.5,0
set style line 1 lc rgb "red"   lw 3.0 lt 2 # red solid
set style line 2 lc rgb "black" lw 2.0 lt 0 # black dashed
set style line 6 lc rgb "black" lw 5.5 lt 0 # black dashed
set style line 3 lc rgb "black"  lw 3.0 lt 1 # blue solid
set style line 4 lc rgb "dark-violet" lw 3.0 lt 4 # black dotted
set style line 5 lc rgb "blue" lw 3.0 lt 4 # black dotted

set label '$Q^2>100$ GeV$^2$' at 0.07,0.21
set label '$15<p_T^\nu<20$ GeV' at 0.07,0.1
set label '$0.1<y<0.85$' at 0.07,0.046

set table 'up_smooth'
plot 'data/zh_theory_new.dat' using 1:3 smooth sbezier

set table 'down_smooth'
plot 'data/zh_theory_new.dat' using 1:4 smooth sbezier
unset table

plot '< paste up_smooth down_smooth' u 1:2:5 with filledcurves fillstyle solid 0.5 lw 4.5 lc rgb 'orange' title '',\
     'data/zh_pythia_pm.dat' w lines lc 'blue' lw 5 title '$\textsc{Pythia8}$',\
     'data/zh_theory_new.dat' u 1:2 w lines smooth sbezier lc 'orange' lw 5 title 'NLL',\
#     'data/zh_theory_new.dat' u 1:3 w lines smooth sbezier lc 'orange' lw 5 title 'theory',\
#     'data/zh_theory_new.dat' u 1:4 w lines smooth sbezier lc 'orange' lw 5 title 'theory',\
#     'zh_1008.dat' u 1:($2/$3) w lines smooth sbezier lc 'orange' lw 2.5 title 'Theory',\
#     'fig5_ave.dat' u 1:2 lc 'red' lw 5 title '',\
#     'y0.dat' w lines lc 'black' title ''
#     './j3_R1_frac.dat' u 1:($2*100) smooth sbezier lt 7 lw 5 title '$0.5 < z_\Lambda < 0.8$',\
#     './0125_zh_5to8.dat' u 1:2 smooth sbezier lt 6 lw 5 title '$0.5 < z_h < 0.8$',\
#     './zh_1to5.dat' u 1:2 w points ls 3 title '',\
#     './zh_5to8.dat' u 1:2 w points ls 4 title '',\
#     './teq_zh_1to5.dat' u 1:2 notitle,\
#     './cteq_zh_5to8.dat' u 1:2 notitle,\

unset ylabel
unset label
set ytics format ""

unset multiplot
set output # finish the current output file
system('latex zh_theory.tex && dvips zh_theory.dvi && ps2eps -f zh_theory.ps && epstopdf zh_theory.eps && rm zh_theory.inc.eps zh_theory.ps *.aux *.dvi *.log *.tex')