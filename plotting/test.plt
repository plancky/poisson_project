set xlabel "x"
set ylabel "y"
set key top
#set border 4095
#set term X11
set xrange [0:4]
set yrange [0:4.4]
#set zrange [-.25:1]
set samples 80,88
set isosamples 80,88
set view 50,50
set title "pm3d demo. Radial sinc function. Default options."
set pm3d at bs; set palette rgbformulae 30,31,32
#show pm3d
#show palette
unset surface
#set view map
#set dgrid3d 30,30 gauss .75,0.5
#set dgrid3d 400,440
splot '/home/planck/Desktop/secc_project_proposal/dats/dataset4_e-8/sor_poisson2d_0.05_1.9_1.dat' 
#splot sin(sqrt(x**2+y**2))/sqrt(x**2+y**2)


pause -1 "Hit any key to continue"
