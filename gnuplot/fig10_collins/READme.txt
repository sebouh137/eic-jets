There are 12 data files and one folder in this folder:

--------------------------------------------------------
For the upper plot in Fig.3 (Electron beam)
Left panel:
-> fig3_pip_e-_p1.dat (orange curve)
-> fig3_pim_e-_p1.dat (blue curve)
Middle panel:
-> fig3_pip_e-_p2.dat (orange curve)
-> fig3_pim_e-_p2.dat (blue curve)
Right panel:
-> fig3_pip_e-_p3.dat (orange curve)
-> fig3_pim_e-_p3.dat (blue curve)

In each file "*.dat":
The first column is zh;
the second column is the central value of A_UT;
the third colunm is the uncertainty value of A_UT.

--------------------------------------------------------
For the lower plot in Fig.3 (Positron beam):
Left panel:
-> fig3_pip_e+_p1.dat (orange curve)
-> fig3_pim_e+_p1.dat (blue curve)
Middle panel:
-> fig3_pip_e+_p2.dat (orange curve)
-> fig3_pim_e+_p2.dat (blue curve)
Right panel:
-> fig3_pip_e+_p3.dat (orange curve)
-> fig3_pim_e+_p3.dat (blue curve)

In each file "*.dat":
The first column is zh;
the second column is the central value of A_UT;
the third colunm is the uncertainty value of A_UT.

--------------------------------------------------------
In the folder "gnuplot/":
-> Tow gnuplot scripts are given for making Fig.3. 
-> To make the plot, run:
$ gnuplot fig3_upper.x 
$ gnuplot fig3_lower.x
