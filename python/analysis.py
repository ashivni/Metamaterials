import numpy
import matplotlib as mpl 
import hierarchicalMaterials as hm
import scipy.optimize 

def stress_conc(nx=100,ny=70,l0=1,levels=0):
	fig = mpl.pyplot.figure()
	fig.subplots_adjust(bottom=0.2,left=0.2)
	ax = fig.add_subplot(111)
	ax.set_xscale('log')
	ax.set_yscale('log')
	"""
	#ax.axis('off')
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	for tic in ax.xaxis.get_major_ticks():
		tic.tick1On = tic.tick2On = False
	for tic in ax.yaxis.get_major_ticks():
		tic.tick1On = tic.tick2On = False
	"""

	ax.set_aspect('equal')

	notch_len = []
	stress_conc = []
	notch_len_frac = numpy.linspace(1.0/8,1.0/2,5)

	for nlf in notch_len_frac:
		print 'Solving ', nlf
		hg = hm.hierarchical_grid(nx,ny,levels,notch=True,l0=l0,notch_len=nlf)
		hg._solve()
		notch_len.append(hg._level_notch_len[0])
		stress_conc.append(hg._level_eqns[0]['curr'].max())

	notch_len = numpy.array(notch_len)
	stress_conc = numpy.array(stress_conc)

	p = scipy.optimize.leastsq(stress_conc_fit_diff,numpy.zeros(2), args = (notch_len, stress_conc))
	C = p[0]
	line, = ax.plot(notch_len,stress_conc,'ko')
	line, = ax.plot(notch_len,stress_conc_fit_func(C,notch_len),'k-')

	mpl.pyplot.draw()

	return numpy.array(notch_len), numpy.array(stress_conc)


def stress_conc_fit_func(C, nl):
	return C[0]*(nl**0.5) + C[1]

def stress_conc_fit_diff(C,nl,sc):
	return sc - stress_conc_fit_func(C,nl)
