set title "Variation of iterations with omega h=0.1, rtol=1e-8"
set xlabel "\omega "
set ylabel "No. of iterations for convergence, rtol = e-8"
set terminal pngcairo
set terminal png font arial 30 size 2400,1300
set xrange [1:2]
set ytics 200
set xtics 0.1
set grid xtics,ytics 
set output 'omega.png'
plot '/home/planck/Desktop/secc_project_proposal/omega_variation.dat' notitle w lp
