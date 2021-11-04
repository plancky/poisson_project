set title "splot with \"set pm3d\" (implemented with some terminals)"
#set pm3d
set hidden3d
set auto
set isosamples 60
#splot '/home/planck/Desktop/secc_project_proposal/dats/dataset1_h2/sor_poisson2d_0.1_1_1.dat' 
#splot '/home/planck/Desktop/secc_project_proposal/dats/dataset2_e-8/sor_poisson2d_0.01_1.8979_1.dat' 
plot '/home/planck/Desktop/secc_project_proposal/dats/dataset2_e-8/sor_poisson2d_0.01_1.8979_1.dat' with image
pause -1 "Hit any key to continue"