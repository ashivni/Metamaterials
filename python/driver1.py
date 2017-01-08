import hierarchicalMaterials as hm
import pylab
import pickle
import numpy

magnification = [8, 16, 24, 32, 40, 48, 56, 64]
nx, ny = 17, 10
max_stress = []
notch_len = []
density = []
toughness = []

for m in magnification:
	print m
	hg1 = hm.hierarchical_grid(nx,ny,1,magnification=m,l0=1,notch=True,damage=False)
	hg1._solve()
	print "Plotting"
	pylab.ioff()
	hg1.plot(show_lines=False,hide_axis=True,show_broken=False,show_nodes=False,show_triangles=True)
	pylab.ion()

	hg1.simulate_fracture()
	max_stress.append(max(hg1._level_stress[0]))
	notch_len.append(hg1._level_notch_len[0])
	density.append(hg1._level_density[0])
	toughness.append(max_stress[-1]*(notch_len[-1]**0.5))

	res = {}
	res['mag'] = numpy.array(magnification)
	res['stress'] = numpy.array(max_stress)
	res['notch_len'] = numpy.array(notch_len)
	res['density'] = numpy.array(density)
	res['toughess'] = numpy.array(toughness)

	f = open('Nx_%d_Ny_%d_NotchData.pic'%(nx,ny),'w+')
	pickle.dump(res,f)
	f.close()

