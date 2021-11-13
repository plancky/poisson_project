set title "Successive Over-relaxation Method h = 0.0125"
set xlabel "x (µm)"
set ylabel "y (µm)"
set zlabel "U (Volts)"
set terminal pngcairo
set terminal png font arial 32 size 2400,1300
set xrange [0:4]
set yrange [0:4.4]
set samples 80,88
set isosamples 80,88
set view 50,50
set pm3d at bs; set palette rgbformulae 30,31,32
unset surface
set view map
set output 'sor_map.png'
splot '/home/planck/Desktop/spp/dats/dataset4_e-8/sor_poisson2d_0.0125_1_1.dat' notitle 
