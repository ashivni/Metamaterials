import hierarchicalMaterials as hm
import pylab
import pickle
import numpy

magnification = [64]
damage = [0.1,0.2,0.3]
nx, ny = 17, 10
nSim = 100
for m in magnification:
	res = {}
	for p in damage:
		max_stress = []
		for i in range(nSim):
			print m, p, i
			hg1 = hm.hierarchical_grid(nx,ny,1,magnification=m,l0=1,notch=False,damage=True,damage_fraction=p)
			hg1.simulate_fracture()
			
			max_stress.append(max(hg1._level_stress[0]))

			res['mag_%d_damage_%.2f_stress'%(m,p)] = numpy.array(max_stress)

			f = open('Nx_%d_Ny_%d_Mag_%d_StrengthData.pic'%(nx,ny,m),'w+')
			pickle.dump(res,f)
			f.close()

