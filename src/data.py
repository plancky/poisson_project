FILE = '/home/planck/Desktop/secc_project_proposal/dats/dataset4_e-8/sor_poisson2d_0.05_1.9_1.dat'
f = open(FILE,'r')
prev_x = f.readline().split()[0]
lst1= ""
for line in f :
    new_x = line.split(' ')[0]
    if new_x!=prev_x:
        lst1+="\n"+line
        prev_x = new_x
    else:
        lst1+=line
f.close()
f2 =open(FILE,"w")
f2.write(lst1)
f2.close()
