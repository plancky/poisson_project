set title "3D surface plot of Potential distribution inside the capacitor"
set xlabel "x (µm)"
set ylabel "y (µm)"
set zlabel "U (Volts) "
set terminal pngcairo
set key top
set xrange [0:4]
set yrange [0:4.4]
#set zrange [-.25:1]
set samples 80,88
set isosamples 80,88
set view 50,50
set pm3d at bs; set palette rgbformulae 30,31,32
set output sprintf('finished.png')
#set dgrid3d 80,88 gauss .75,0.5
unset surface
#set contour 
#set cntrparam levels 15
#set view map


list = system('ls ../dats/dataset4_e-8/')

do for [file in list] {
    set output sprintf('%s.png', file)
    set title sprintf("%s", file)
    splot sprintf("../dats/dataset4_e-8/%s",file)
}
#splot '/home/planck/Desktop/secc_project_proposal/dats/dataset4_e-8/sor_poisson2d_0.05_1.9_1.dat' 
pause -1 "Hit any key to continue"
