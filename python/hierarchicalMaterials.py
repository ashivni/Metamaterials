import numpy
import pickle
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg
import matplotlib as mpl 
import matplotlib.cm as cm
import collections
import sys

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

class node:
	def __init__(self, i, j,_all_nodes=None,pad=9):
		""" Node at lattice index i, j. The basis vectors are v_1 = e_1/2, v_2 = sqrt(3)*0.5 e_2 """
		fmt_str = '_%0' + str(pad) + 'd_%0' + str(pad) + 'd'
		self._id = fmt_str%(i,j)
		if _all_nodes is not None and self._id in _all_nodes:
			raise ValueError('Node at index (%d,%d) already exists'%(i,j))
		elif _all_nodes is not None:
			_all_nodes.add(self._id)
		self._i = i
		self._j = j
		self._bonded_to = []

	def bond_to(self, node):
		if node._id in self._bonded_to and self._id in node._bonded_to:
			raise ValueError("[Attempting Bond] Nodes %s, %s already bonded."%(self._id, node._id))
		elif node._id in self._bonded_to or self._id in node._bonded_to:
			raise ValueError("[Attempting Bond] Nodes %s, %s partially bonded."%(self._id, node._id))
		else:
			self._bonded_to.append(node._id)
			node._bonded_to.append(self._id)

	def is_bonded_to(self,node):
		return node._id in self._bonded_to

	def debond_from(self,node):
		if node._id not in self._bonded_to and self._id not in node._bonded_to:
			raise ValueError("[Attempting debond] Nodes %s, %s not bonded."%(self._id, node._id))
		elif node._id not in self._bonded_to or self._id not in node._bonded_to:
			raise ValueError("[Attempting debond] Nodes %s, %s partially bonded."%(self._id, node._id))
		else:
			#print "Removing ", node._id, " from ", str(self._bonded_to)
			self._bonded_to.remove(node._id)
			#print "Removed: ", str(self._bonded_to)
			
			#print "Removing ", self._id, " from ", str(node._bonded_to)
			node._bonded_to.remove(self._id)
			#print "Removed: ", str(node._bonded_to)

	def coordinates(self):
		return 0.5*self._i, 0.5*self._j*(3.0**0.5)

	def global_node_number(self,nx):
		#FIXME
		return (self._j*nx + self._i)*2

	def is_connected_to(self, target, node_dict,verbose=False):
		""" 
		Does a BF-Search (starting at self) to check connectivity with the given node 
		Worse case runtime is order(N), though for nearby connected nodes should typically run 
		much faster
		"""


		dq = collections.deque()
		if verbose:
			print "Starting search at ", self._id, ". Start neigh list = ", str(self._bonded_to)
			print "Q Len ", len(dq)
		dq.append(self)
		visited = set()
		queued = set()

		while len(dq) > 0:
			n = dq.popleft()
			if n is target:
				return True
			else:
				visited.add(n._id)

			for ng_id in n._bonded_to:
				if ng_id not in visited and ng_id not in queued:
					dq.append( node_dict[ng_id])
					queued.add(ng_id)

		return False

	def __str__(self):
		_c = self.coordinates()
		return 'Node [ ' + '%d, %d, %.2f, %.2f '%(self._i, self._j, _c[0], _c[1]) + ' ]' 
	
	def _node2str(self):
		_c = self.coordinates()
		return ['%d, %d, %.2f, %.2f'%(self._i, self._j,_c[0],_c[1])]

class bond:
	def __init__(self, n1, n2, bond_type,_all_bonds=None):
		""" Bond between nodes n1, n2 """
		self._id = n1._id + '_' + n2._id
		if _all_bonds is not None and self._id in _all_bonds:
			raise ValueError('Bond between nodes %s, %s already exists'%(n1._id,n2._id))
		elif _all_bonds is not None:
			_all_bonds.add(self._id)

		self._n1 = n1
		self._n2 = n2
		self._n1.bond_to(self._n2)
		self._bond_type = bond_type

	def coordinates(self,nx,ny):
		""" returns tupples of plot-able line segments corresponding to the bond """
		_c1, _c2 = self._n1.coordinates(), self._n2.coordinates()
		rt3 = 3.0**0.5
		if self._bond_type in ['R', 'U', 'UR','I']:
			_X = [(_c1[0], _c2[0])]
			_Y = [(_c1[1], _c2[1])]
		elif self._bond_type == 'PR':
			_X = [(_c1[0], (_c1[0]+nx+_c2[0])*0.5),			( (_c1[0] -nx+_c2[0])*0.5, _c2[0])]
			_Y = [(_c1[1], (_c1[1]+_c2[1])*0.5),			( (_c1[1]+_c2[1])*0.5, _c2[1])]

		elif self._bond_type == 'PU':
			_X = [(_c1[0], (_c1[0]+_c2[0])*0.5),			( (_c1[0]+_c2[0])*0.5, _c2[0])]
			_Y = [(_c1[1], (_c1[1]+_c2[1]+ny*rt3)*0.5),			( (_c1[1]+_c2[1]-ny*rt3)*0.5, _c2[1])]
			
		elif self._bond_type == 'C':
			_X = [(_c1[0], (_c1[0]+nx+_c2[0])*0.5),			( (_c1[0] -nx+_c2[0])*0.5, _c2[0])]
			_Y = [(_c1[1], (_c1[1]+_c2[1]+ny*rt3)*0.5),			( (_c1[1]+_c2[1]-ny*rt3)*0.5, _c2[1])]
		
		elif self._bond_type == 'PUR':
			_X = [(_c1[0], (_c1[0]+nx+_c2[0])*0.5),			( (_c1[0] -nx+_c2[0])*0.5, _c2[0])]
			_Y = [(_c1[1], (_c1[1]+_c2[1])*0.5),			( (_c1[1]+_c2[1])*0.5, _c2[1])]
		
		return _X, _Y


	def refine(self,magnification):
		_c1, _c2 = self._n1.coordinates(), self._n2.coordinates()
		_t_vec = [ (self._n2._i - self._n1._i)/magnification, (self._n2._j - self._n1._j)/magnification]
		_t_vec_real = [(_c2[0] - _c1[0])/magnification, (_c2[1] - _c1[1])/magnification]	# Basic translation vector
		_l = int(round((_t_vec_real[0]**2.0 + _t_vec_real[1]**2.0)**0.5))
		_nodes = {}
		_bonds = {}
		for k in range(magnification):
			h = hexagon(self._n1._i + k*_t_vec[0], self._n1._j + k*_t_vec[1],_l)
			for n in h._nodes:
				if n._id not in _nodes:
					_nodes[n._id] = n

			for b in h._bonds:
				if b._id not in _bonds:
					_bonds[b._id] = b
					try:
						_nodes[b._n1._id].bond_to(_nodes[b._n2._id])
					except ValueError:
						pass

		"""
		fig = mpl.pyplot.figure()
		fig.subplots_adjust(bottom=0.2,left=0.2)
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')
		#ax.axis('off')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		for tic in ax.xaxis.get_major_ticks():
			tic.tick1On = tic.tick2On = False
		for tic in ax.yaxis.get_major_ticks():
			tic.tick1On = tic.tick2On = False
		
		ax.set_aspect('equal')

		x, y = [], []
		for n_id, n in _nodes.iteritems():
			_c = n.coordinates()
			x.append(_c[0])
			y.append(_c[1])
		line, = ax.plot(x,y,'ko')

		for b_id, b in _bonds.iteritems():
			_tX,_tY = b.coordinates(_l*10,_l*10)
			for _X, _Y in zip(_tX,_tY):
				line, = ax.plot(_X,_Y,'b-')
		
		mpl.pyplot.draw()
		"""

		return _nodes, _bonds

	def sever(self):
		self._n1.debond_from(self._n2)

	def __str__(self):
		return 'Bond: ' + str(self._n1) + ' <---> ' + str(self._n2)

	def _bond2str(self):
		_sn1 = self._n1._node2str()
		_sn2 = self._n2._node2str()

		return [_sn1[0] + ' <---> ' + _sn2[0] + ' ' + self._bond_type]

	def length(self,nx,ny):
		_c1, _c2 = self._n1.coordinates(), self._n2.coordinates()
		rt3 = 3.0**0.5
		if self._bond_type in ['R', 'U', 'UR','I']:
			return ( (_c1[0] - _c2[0] - nx)**2.0 + (_c1[1] - _c2[1])**2.0)**0.5
		
		elif self._bond_type == 'PR':
			return ( (_c1[0] - _c2[0])**2.0 + (_c1[1] - _c2[1] - ny*rt3 )**2.0)**0.5

		elif self._bond_type == 'PU':
			return ( (_c1[0] - _c2[0])**2.0 + (_c1[1] - _c2[1] - ny*rt3 )**2.0)**0.5
			
		elif self._bond_type == 'PUR':
			return ( (_c1[0] - _c2[0]+nx)**2.0 + (_c1[1] - _c2[1] - ny*rt3 )**2.0)**0.5
	
class hexagonal_unit_cell:
	def __init__(self, a, i, j,_all_cells=None):
		self._id = '%d_%d'%(i,j)
		if _all_cells is not None and self._id in _all_cells:
			raise ValueError("Unit cell at (%d, %d) already exits"%(i,j))
		if _all_cells is not None:	
			_all_cells.add(self._id)

		self._a = a
		self._i = i
		self._j = j

		self._nodes = []
		self._interior_bonds = []

		self._nodes.append(node(self._i*2*self._a, self._j*2*self._a))
		self._nodes.append(node((self._i*2+1)*self._a, (self._j*2+1)*self._a))
		self._interior_bonds.append(bond(self._nodes[0], self._nodes[1],'I'))

		self._neighbors = []
		self._exterior_bonds = []

	def add_neighbor(self, neighbor, neigh_type):

		if neighbor._id in self._neighbors and self._id in neighbor._neighbors:
			raise ValueError("[Attempting neighboring] Unit cells %s and %s are already neighbors"%(self._id, neighbor._id))
		elif neighbor._id in self._neighbors or self._id in neighbor._neighbors:
			raise ValueError("[Attempting neighboring] Unit cells %s and %s are partial neighbors"%(self._id, neighbor._id))
		else:
			if neigh_type == 'R' or neigh_type == 'PR':
				bonds = [(0,0), (1,1), (1,0)]
			elif neigh_type == 'U' or neigh_type == 'PU':
				bonds = [(1,0)]
			elif neigh_type == 'UR' or neigh_type == 'PUR':
				bonds = [(1,0)]
			elif neigh_type == 'C':
				bonds = [(1,0)]
			else:
				raise ValueError("Unknown neighboring cell type")

			for b in bonds:
				_b = bond(self._nodes[b[0]],neighbor._nodes[b[1]], neigh_type)
				self._exterior_bonds.append(_b)
				neighbor._exterior_bonds.append(_b)
				self._neighbors.append(neighbor._id)
				neighbor._neighbors.append(self._id)

	def _nodes2str(self):
		_s = []
		for n in self._nodes:
			_s.extend(n._node2str())

		return _s

	def _bonds2str(self):
		_s = []
		for b in self._interior_bonds:
			_s.extend(b._bond2str())
		
		for b in self._exterior_bonds:
			_s.extend(b._bond2str())

		return _s

class hexagon:
	def __init__(self, i, j, a):
		""" Hexagon of side 'a' centered at i, j """
		self._i = i
		self._j = j
		self._a = a
		self._nodes = []
		self._bonds = []
		self._nodes.append(node(self._i,self._j))
		self._nodes.append(node(self._i+2*self._a, self._j))
		self._nodes.append(node(self._i+self._a, self._j+self._a))
		self._nodes.append(node(self._i-self._a, self._j+self._a))
		self._nodes.append(node(self._i-2*self._a, self._j))
		self._nodes.append(node(self._i-self._a, self._j-self._a))
		self._nodes.append(node(self._i+self._a, self._j-self._a))
		for k in range(6):
			self._bonds.append(bond(self._nodes[0],self._nodes[k+1],'I'))

		for k in range(6):
			self._bonds.append(bond(self._nodes[1+(k+1)%6],self._nodes[1+(k+2)%6],'I'))


	def node_coordinates(self):
		_X, _Y = [], []
		for node in self._nodes:
			_x, _y = node.coordinates()
			_X.append(_x)
			_Y.append(_y)

		return _X, _Y

	def plot(self,lines=False):
		fig = mpl.pyplot.figure()
		fig.subplots_adjust(bottom=0.2,left=0.2)
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')
		#ax.axis('off')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		for tic in ax.xaxis.get_major_ticks():
			tic.tick1On = tic.tick2On = False
		for tic in ax.yaxis.get_major_ticks():
			tic.tick1On = tic.tick2On = False
		
		ax.set_aspect('equal')

		x, y = self.node_coordinates()
		line, = ax.plot(x,y,'ko')

		if lines:
			for _b in self._bonds:
				_tX,_tY = _b.coordinates(self._a*self._magnification,self._a*self._magnification)
				for _X, _Y in zip(_tX,_tY):
					line, = ax.plot(_X,_Y,'b-')
		
		mpl.pyplot.draw()

class hexagonal_grid:
	def __init__(self, nx, ny, level,periodic=False,**kwargs):
		self._nx = nx
		self._ny = ny
		self._level = level
		self._l0 = kwargs.get('l0',1)
		self._magnification = kwargs.get('magnification',10)
		self._a = self._l0*(self._magnification**self._level)
		self._master_grid = {}
		self._all_cells = set()
		for i in xrange(self._nx):
			for j in xrange(self._ny):
				cell_id = (j*self._nx + i)*self._a
				self._master_grid[cell_id] = hexagonal_unit_cell(self._a, i, j, self._all_cells)
	
		# Add neighbors
		for i in xrange(self._nx):
			for j in xrange(self._ny):
				cell_id = (i + j*self._nx)*self._a
				# offsets in i, j for right, upper, and upper-right neighbors
				for o_i, o_j, neigh_type in zip([1,0,1],[0,1,1],['R','U','UR']):
					if o_i + i == self._nx and o_j + j == self._ny:
						neigh_type = 'C'
					elif o_i + i == self._nx and o_j == 1:
						neigh_type = 'PUR'
					elif o_i + i == self._nx:
						neigh_type = 'PR'
					elif o_j + j == self._ny:
						neigh_type = 'PU'
					n_i, n_j= (i + o_i)%self._nx, (j + o_j)%self._ny
					neigh_id = (n_i + n_j*self._nx)*self._a
					if periodic or not periodic and neigh_type in ['R','U','UR']:
						self._master_grid[cell_id].add_neighbor(self._master_grid[neigh_id], neigh_type)

	def node_dict(self):
		_n_dict = {}
		for cell_id, cell in self._master_grid.iteritems():
			for _n in cell._nodes:
				if _n._id in _n_dict:
					raise ValueError("Duplicate node")
				else:
					_n_dict[_n._id] = _n
		
		return _n_dict

	def bond_dict(self):
		_b_dict = {}
		for cell_id, cell in self._master_grid.iteritems():
			for _b in cell._interior_bonds:
				if _b._id in _b_dict:
					pass
				else:
					_b_dict[_b._id] = _b
			
			for _b in cell._exterior_bonds:
				if _b._id in _b_dict:
					pass
				else:
					_b_dict[_b._id] = _b

		return _b_dict

	def node_coordinates(self):
		_X, _Y = [], []
		for cell_id, cell in self._master_grid.iteritems():
			for node in cell._nodes:
				_x, _y = node.coordinates()
				_X.append(_x)
				_Y.append(_y)

		return _X, _Y

	def _nodes2str(self,coords=True):
		_s = []
		for cell_id, cell in self._master_grid.iteritems():
			_s.extend(cell._nodes2str())

		_s.sort()
		return _s

	def _bonds2str(self):
		_s = []
		for cell_id, cell in self._master_grid.iteritems():
			_s.extend(cell._bonds2str())
		
		_s.sort()

		return _s
		

	def plot(self,lines=False):
		fig = mpl.pyplot.figure()
		fig.subplots_adjust(bottom=0.2,left=0.2)
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')
		#ax.axis('off')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		for tic in ax.xaxis.get_major_ticks():
			tic.tick1On = tic.tick2On = False
		for tic in ax.yaxis.get_major_ticks():
			tic.tick1On = tic.tick2On = False
		
		ax.set_aspect('equal')
	
		x, y = self.node_coordinates()
		line, = ax.plot(x,y,'k.')


		if lines:
			for cell_id, cell in self._master_grid.iteritems():
				for _ib in cell._interior_bonds:
					_tX,_tY = _ib.coordinates(self._nx*self._a,self._ny*self._a)
					for _X, _Y in zip(_tX,_tY):
						line, = ax.plot(_X,_Y,'r-')
			
				for _eb in cell._exterior_bonds:
					_tX,_tY = _eb.coordinates(self._nx*self._a,self._ny*self._a)
					for _X, _Y in zip(_tX,_tY):
						if _eb._bond_type in ['U','R','UR']:
							line, = ax.plot(_X,_Y,'b-')
						elif _eb._bond_type in ['PU','PR','PUR','C']:
							line, = ax.plot(_X,_Y,'b--')
		
		mpl.pyplot.draw()
	 
		
class hierarchical_grid:
	def __init__(self,nx,ny,levels,**kwargs):
		self._n_step = 0
		self._nx = nx
		self._ny = ny
		self._levels = levels
		self._l0 = kwargs.get('l0',1)
		self._magnification = kwargs.get('magnification',10)
		self._a = self._l0*(self._magnification**self._levels)
		self._ly = self._a*(self._ny*2 -1)*(3.0**0.5)*0.5
		self._lx = self._a*(self._nx*2 -1)*0.5
		self._ly_ind = self._a*(self._ny*2 -1)
		self._lx_ind = self._a*(self._nx*2 -1)
		self._notch = kwargs.get('notch',False)
		self._damage = kwargs.get('damage',False)
		self._damage_fraction = kwargs.get('damage_fraction',0.0)
		self._outline = hexagonal_grid(self._nx,self._ny,self._levels,periodic=False,l0=self._l0,magnification=self._magnification)
		self._level_nodes, self._level_bonds, self._level_broken_bonds = {}, {}, {}
		self._level_not_broken, self._level_is_solved, self._level_is_build = {}, {}, {}
		self._level_voltage, self._level_current = {}, {}
		self._level_strain, self._level_strain_evp, self._level_scaled_stress, self._level_stress = {}, {}, {}, {}
		self._level_density = {}
		self._level_nodes[self._levels] = self._outline.node_dict()
		self._level_bonds[self._levels] = self._outline.bond_dict()
		self._level_broken_bonds[self._levels] = []
		
		self._level_not_broken[self._levels] = True
		self._level_is_build[self._levels] = False
		self._level_is_solved[self._levels] = False
		self._level_voltage[self._levels] = []
		self._level_current[self._levels] = []
		self._level_strain[self._levels] = []
		self._level_strain_evp[self._levels] = []
		self._level_scaled_stress[self._levels] = []
		self._level_stress[self._levels] = []
		self._level_density[self._levels] = len(self._level_bonds[self._levels])/(self._lx*self._ly)
		_l = self._levels -1
		while _l >= 0:
			_nodes, _bonds = {}, {}

			for _b_up_id, _b_up in self._level_bonds[_l+1].iteritems():

				_n_refined, _b_refined = _b_up.refine(magnification=self._magnification)

				for _n_id, _n in _n_refined.iteritems():
					if _n_id not in _nodes:
						_nodes[_n_id] = _n
				
				for _b_id, _b in _b_refined.iteritems():
					if _b_id not in _bonds:
						_bonds[_b_id] = _b
					
					_n1 = _nodes[_b._n1._id]
					_n2 = _nodes[_b._n2._id]

					if not _n1.is_bonded_to(_n2):
						_n1._bonded_to.append(_n2._id)

					if not _n2.is_bonded_to(_n1):
						_n2._bonded_to.append(_n1._id)


			self._level_nodes[_l] = _nodes
			self._level_bonds[_l] = _bonds

				
			for _bond_id, _bond in self._level_bonds[_l].iteritems():
				_n1, _n2 = _bond._n1, _bond._n2
				_bond._n1 , _bond._n2 = self._level_nodes[_l][_n1._id], self._level_nodes[_l][_n2._id]

			_duplicate_bonds = set()
			for _bond_id, _bond in self._level_bonds[_l].iteritems():
				_n1, _n2 = _bond._n1, _bond._n2
				if _bond_id not in _duplicate_bonds:
					_dup_id = _n2._id + '_' + _n1._id
					if _dup_id in self._level_bonds[_l]:
						_duplicate_bonds.add(_dup_id)

			for _dup_id in _duplicate_bonds:
				self._level_bonds[_l].pop(_dup_id)
			
					
			self._level_not_broken[_l] = True
			self._level_broken_bonds[_l] = []
			self._level_is_solved[_l] = False
			self._level_is_build[_l] = False 
			self._level_voltage[_l] = []
			self._level_current[_l] = []
			self._level_strain[_l] = []
			self._level_strain_evp[_l] = []
			self._level_scaled_stress[_l] = []
			self._level_stress[_l] = []
			self._level_density[_l] = len(self._level_bonds[_l])/(self._lx*self._ly)

			_l -= 1


		if self._notch:
			_l = self._levels
			self._level_notched_nodes = {}
			self._level_notched_bonds = {}
			self._level_notch_len = {}
			while _l >= 0:


				_i_min, _j_min = numpy.inf, numpy.inf
				_i_max, _j_max = -numpy.inf, -numpy.inf
				

				for node_id in self._level_nodes[_l]:
					x, y = [int(z) for z in node_id.split('_')[1:]]
					_i_min = _i_min if _i_min < x else x
					_i_max = _i_max if _i_max > x else x
					_j_min = _j_min if _j_min < y else y
					_j_max = _j_max if _j_max > y else y

				_notch_len_ind = kwargs.get('notch_len',0.25)*(_i_max - _i_min)
				_notch_len_ind = numpy.round(_notch_len_ind) + self._l0*(self._magnification**self._levels)
				
				_j_mid = 0.5*(_j_min + _j_max)
				if abs( round(_j_mid) - _j_mid)  < 1E-3:
					_j_mid += 0.25

				_notched_nodes = {}
				_notched_bonds = {}

				for node_id in self._level_nodes[_l]:
					node_x = int(node_id.split('_')[1])
					node_y = int(node_id.split('_')[2])
					neighs = self._level_nodes[_l][node_id]._bonded_to
					for neigh_id in neighs:
						neigh_x = int(neigh_id.split('_')[1])
						neigh_y = int(neigh_id.split('_')[2])
						if node_x - _i_min <= _notch_len_ind and (node_y - _j_mid)*(neigh_y - _j_mid) < 0:
							self._level_nodes[_l][node_id].debond_from(self._level_nodes[_l][neigh_id])
							if node_id + '_' + neigh_id in self._level_bonds[_l]:
								_notched_bonds[node_id + '_' + neigh_id] = self._level_bonds[_l].pop(node_id + '_' + neigh_id)
							if neigh_id + '_' + node_id in self._level_bonds[_l]:
								_notched_bonds[neigh_id + '_' + node_id] = self._level_bonds[_l].pop(neigh_id + '_' + node_id)

							if node_id not in _notched_nodes:
								_notched_nodes[node_id] = set()
							_notched_nodes[node_id].add(neigh_id)
							
							if neigh_id not in _notched_nodes:
								_notched_nodes[neigh_id] = set()
							_notched_nodes[neigh_id].add(node_id)
				
				self._level_notched_nodes[_l] = _notched_nodes
				self._level_notched_bonds[_l] = _notched_bonds
				self._level_notch_len[_l] = _notch_len_ind*0.5

				_l -= 1
		
		if self._damage:
			_l = self._levels

			self._level_damaged_nodes = {}
			self._level_damaged_bonds = {}
			while _l >= 0:
				# Test neigh lists
				for _node_id, _node in self._level_nodes[_l].iteritems():
					for _neigh_id in _node._bonded_to:
						if _node_id not in self._level_nodes[_l][_neigh_id]._bonded_to:
							raise Exception('Bad neighbor list')
				
				# Test for object consistency
				for _bond_id, _bond in self._level_bonds[_l].iteritems():
					_n1, _n2 = _bond._n1, _bond._n2
					_nl1 , _nl2 = self._level_nodes[_l][_n1._id], self._level_nodes[_l][_n2._id]
					if _n1 is not _nl1 or _n2 is not _nl2:
						raise Exception('Dual identities')

				_damaged_nodes = {}
				_damaged_bonds = {}
				

				for _bond_id, _bond in self._level_bonds[_l].iteritems():
					_n1, _n2 = _bond._n1, _bond._n2
					node_x = int(_n1._id.split('_')[1])
					node_y = int(_n1._id.split('_')[2])
					if numpy.random.rand() < self._damage_fraction and _n1._id + '_' + _n2._id not in _damaged_bonds and _n2._id + '_' + _n1._id not in _damaged_bonds:
						_n1.debond_from(_n2)
						if _n1.is_connected_to(_n2, self._level_nodes[_l]):

							if _n1._id not in _damaged_nodes:
								_damaged_nodes[_n1._id] = set()
							_damaged_nodes[_n1._id].add(_n2._id)
							
							if _n2._id not in _damaged_nodes:
								_damaged_nodes[_n2._id] = set()
							_damaged_nodes[_n2._id].add(_n1._id)

							_damaged_bonds[_bond_id] = _bond
						else:
							_n1.bond_to(_n2)

					
				for _damaged_bond_id in _damaged_bonds:
					self._level_bonds[_l].pop(_damaged_bond_id)

				self._level_damaged_nodes[_l] = _damaged_nodes
				self._level_damaged_bonds[_l] = _damaged_bonds
				for _node_id, _node in self._level_nodes[_l].iteritems():
					for _neigh_id in _node._bonded_to:
						if _node_id not in self._level_nodes[_l][_neigh_id]._bonded_to:
							print "Node ", _node_id, " Neigh list: ", str(_node._bonded_to)
							print "Node ", _neigh_id, " Neigh list: ", str(self._level_nodes[_l][_neigh_id]._bonded_to)
							raise Exception('Bad neighbor list ', _node_id, _neigh_id)
				
				for _node_id, _node in self._level_nodes[_l].iteritems():
					if len(_node._bonded_to) == 0:
						raise Exception('Isolated Node ' + _node_id)
			
				
				for _damaged_bond_id in _damaged_bonds:
					if _damaged_bond_id in self._level_bonds[_l]:
						raise Exception('Damaged and undamaged bond')

				_l -= 1


		self._EY = kwargs.get('EY',1.0)
		self._level_eqns = {}

	def _solve(self):
		self._build_eqns()
		_l = self._levels 
		EY = self._EY

		while _l >= 0:
			if not self._level_is_solved[_l]:
				_eqns = self._level_eqns[_l]
				L, B, C, V = _eqns['laplacian'], _eqns['load'], _eqns['conductance'], _eqns['vol_diag']
				_N_interior, _N_exterior_up, _N_exterior_dn = _eqns['_N_interior'], _eqns['_N_exterior_up'], _eqns['_N_exterior_dn']

				_int_vol = scipy.sparse.linalg.spsolve(L,B)
				V.data[:_N_interior] = _int_vol

				curr = (V*C - C*V)	# matrix of current between nodes i, j
				curr_up = curr[_N_interior:_N_interior+_N_exterior_dn,:].sum()
				curr_dn = curr[_N_interior+_N_exterior_dn:_N_interior+_N_exterior_dn+_N_exterior_up,:].sum()

				self._level_eqns[_l]['curr'] = curr
				self._level_eqns[_l]['curr_up'] = curr_up
				self._level_eqns[_l]['curr_dn'] = curr_dn
			
				if abs(curr_up) < 1E-2:
					self._level_not_broken[_l] = False
					
				self._level_is_solved[_l] = True
			_l -= 1



	def _check(self,verbose=False):
		self._solve()
		_l = self._levels 

		while _l >= 0:
			if verbose:
				print 'Checking level ', _l

			_eqns = self._level_eqns[_l]
			curr = _eqns['curr']
			V = _eqns['vol_diag']
			_N_interior, _N_exterior_up, _N_exterior_dn = _eqns['_N_interior'], _eqns['_N_exterior_up'], _eqns['_N_exterior_dn']
			_interior_node_ids, _exterior_up_node_ids, _exterior_dn_node_ids = _eqns['interior_node_ids'], _eqns['exterior_up_node_ids'], _eqns['exterior_dn_node_ids']
			_interior_node_indices, _exterior_up_node_indices, _exterior_dn_node_indices = _eqns['interior_node_indices'], _eqns['exterior_up_node_indices'], _eqns['exterior_dn_node_indices']
			N = _N_interior + _N_exterior_dn + _N_exterior_up
		
			if verbose:

				if self._notch:
					print 'Notched Bonds: '
					for _b_id, _b in self._level_notched_bonds[_l].iteritems():
						print '\t (%d,%d) notched from (%d,%d)'%(_b._n1._i, _b._n1._j, _b._n2._i, _b._n2._j)

				print 'Node Voltages: '
				for i in range(N):
					if i < _N_interior:
						node_id = _interior_node_ids[i]
					elif i < _N_interior + _N_exterior_dn:
						node_id = _exterior_dn_node_ids[i-_N_interior]
					else:
						node_id = _exterior_up_node_ids[i-_N_interior-_N_exterior_dn]
					
					node_vol = V[i,i]

					x, y = [int(z) for z in node_id.split('_')[1:]]
					
					print 'Node (%d, %d). Node voltage = %.2f'%(x,y, node_vol)


			up_curr = 0.0
			dn_curr = 0.0
			for i in range(N):

				if i < _N_interior:
					node_id = _interior_node_ids[i]
				elif i < _N_interior + _N_exterior_dn:
					node_id = _exterior_dn_node_ids[i-_N_interior]
				else:
					node_id = _exterior_up_node_ids[i-_N_interior-_N_exterior_dn]
				
				node_vol = V[i,i]

				x, y = [int(z) for z in node_id.split('_')[1:]]
				
				if verbose:
					print 'Checking node (%d, %d, %s). Node voltage = %.2f'%(x,y, node_id, node_vol)

				cur_sum = 0.0
				if verbose and self._notch:
					if node_id in self._level_notched_nodes[_l]:
						for neigh_id in self._level_notched_nodes[_l][node_id]:
							neigh_x, neigh_y = [int(z) for z in neigh_id.split('_')[1:]]
							print bcolors.FAIL + '\t Notched from: (%d,%d)'%(neigh_x, neigh_y) + bcolors.ENDC
							
				for neigh_id in self._level_nodes[_l][node_id]._bonded_to:
					if neigh_id in _interior_node_ids:
						j = _interior_node_indices[neigh_id]
					elif neigh_id in _exterior_dn_node_ids:
						j = _exterior_dn_node_indices[neigh_id] + _N_interior
					elif neigh_id in _exterior_up_node_ids:
						j = _exterior_up_node_indices[neigh_id] + _N_interior + _N_exterior_dn
					else:
						raise ValueError('Unknown neigh_id')

					neigh_vol = V[j,j]
					neigh_curr = curr[i,j]
					cur_sum += neigh_curr
					if verbose:
						neigh_x, neigh_y = [int(z) for z in neigh_id.split('_')[1:]]
						print '\t Neighbor: (%d,%d), vol: %.2f, vol_diff: %.2f, current: %.2f'%(neigh_x, neigh_y, neigh_vol, node_vol - neigh_vol, neigh_curr)

					if abs(neigh_curr - (node_vol - neigh_vol)) > 1E-3:
						print 'Checking node ' + node_id
						print 'Node Voltage = %.2f'%(node_vol)
						print '\t Neigh: %s at voltage %.2f. Vol_diff = %.2f, curr = %.2f'%(neigh_id,neigh_vol,node_vol-neigh_vol, neigh_curr)
						raise Exception('Error')
				
				if verbose:
					print '\t Net current: %.2f'%(cur_sum)

				if i < _N_interior:
					if abs(cur_sum) > 1E-3:
						raise Exception('Nonzero current in interior node')
				elif i < _N_interior + _N_exterior_dn:
					dn_curr += cur_sum
				else:
					up_curr += cur_sum


			if abs(dn_curr + up_curr) > 1E-3:
				raise Exception('Up down mismatch')
			_l -= 1

		print 'All checks passed'
		return True
	
	
	def step(self,**kwargs):
		"""
		Advance simulation one step by breaking one bond at each level. 
		"""

		self._solve()
		_l = self._levels 
		a = self._a
		ly = self._ly
		lx = self._lx
		EY = self._EY

		while _l >= 0:
			if self._level_not_broken[_l]:
				_eqns = self._level_eqns[_l]
				curr = _eqns['curr']
				_N_interior, _N_exterior_up, _N_exterior_dn = _eqns['_N_interior'], _eqns['_N_exterior_up'], _eqns['_N_exterior_dn']
				_interior_node_ids, _exterior_up_node_ids, _exterior_dn_node_ids = _eqns['interior_node_ids'], _eqns['exterior_up_node_ids'], _eqns['exterior_dn_node_ids']
				_interior_node_indices, _exterior_up_node_indices, _exterior_dn_node_indices = _eqns['interior_node_indices'], _eqns['exterior_up_node_indices'], _eqns['exterior_dn_node_indices']
				N = _N_interior + _N_exterior_dn + _N_exterior_up
			
				L = _eqns['laplacian']
				B = _eqns['load']
				C = _eqns['conductance']
				Ci = _eqns['conductance_i']
				Cj = _eqns['conductance_j']

				"""
				# Find which bond to break
				_max_curr, _c_i_break, _c_j_break = csc_argmax(curr)
				print _max_curr, _c_i_break, _c_j_break
				#_max_curr, _c_i_break, _c_j_break = csc_argmax(curr,Ci,Cj)
				print _max_curr, _c_i_break, _c_j_break
				print '***'
				"""
				_max_curr = -numpy.inf
				_n1_id, _n2_id = 0, 0
				for bond_id, bond in self._level_bonds[_l].iteritems():
					n1_b, n2_b = bond._n1._id, bond._n2._id
					if n1_b in _interior_node_indices:
						n1_ind = _interior_node_indices[n1_b]
					elif n1_b in _exterior_dn_node_indices:
						n1_ind = _exterior_dn_node_indices[n1_b] + _N_interior
					else:
						n1_ind = _exterior_up_node_indices[n1_b] + _N_interior + _N_exterior_dn

					if n2_b in _interior_node_indices:
						n2_ind = _interior_node_indices[n2_b] 
					elif n2_b in _exterior_dn_node_indices:
						n2_ind = _exterior_dn_node_indices[n2_b] + _N_interior
					else:
						n2_ind = _exterior_up_node_indices[n2_b] + _N_interior + _N_exterior_dn

					if abs(curr[n1_ind,n2_ind]) > _max_curr:
						_max_curr = abs(curr[n1_ind, n2_ind])
						_n1_id, _n2_id = n1_b, n2_b

				if abs(_max_curr - abs(curr.data).max()) > 1E-3:
					raise Exception('This')
				"""
				if _c_i_break < _N_interior:
					_n1_id = _interior_node_ids[_c_i_break]
				elif _c_i_break < _N_interior + _N_exterior_dn:
					_n1_id = _exterior_dn_node_ids[_c_i_break-_N_interior]
				else:
					_n1_id = _exterior_up_node_ids[_c_i_break-_N_interior-_N_exterior_dn]
				
				if _c_j_break < _N_interior:
					_n2_id = _interior_node_ids[_c_j_break]
				elif _c_j_break < _N_interior + _N_exterior_dn:
					_n2_id = _exterior_dn_node_ids[_c_j_break-_N_interior]
				else:
					_n2_id = _exterior_up_node_ids[_c_j_break-_N_interior-_N_exterior_dn]
				"""
				# Debond
				self._level_nodes[_l][_n1_id].debond_from(self._level_nodes[_l][_n2_id])

				if _n1_id + '_' + _n2_id in self._level_bonds[_l]:
					self._level_broken_bonds[_l].append( self._level_bonds[_l].pop( _n1_id + '_' + _n2_id))
				if _n2_id + '_' + _n1_id in self._level_bonds[_l]:
					self._level_broken_bonds[_l].append( self._level_bonds[_l].pop( _n2_id + '_' + _n1_id))

				# Update matrices 
				_break_node_id = [_n1_id, _n2_id]
				_break_neigh_id = [_n2_id, _n1_id]

				for node_id, neigh_id in zip(_break_node_id,_break_neigh_id):
						
					if neigh_id in _interior_node_indices:
						neigh_ind = _interior_node_indices[neigh_id] 
					elif neigh_id in _exterior_dn_node_indices:
						neigh_ind = _exterior_dn_node_indices[neigh_id] + _N_interior
					elif neigh_id in _exterior_up_node_indices:
						neigh_ind = _exterior_up_node_indices[neigh_id] + _N_interior + _N_exterior_dn
					
					if node_id in _interior_node_indices:
						node_ind = _interior_node_indices[node_id] 
					elif node_id in _exterior_dn_node_indices:
						node_ind = _exterior_dn_node_indices[node_id] + _N_interior
					elif node_id in _exterior_up_node_indices:
						node_ind = _exterior_up_node_indices[node_id] + _N_interior + _N_exterior_dn

					# Update laplacian and load
					if node_id in _interior_node_indices:
						node_x = int(node_id.split('_')[1])
						node_y = int(node_id.split('_')[2])
						diag, load = 0, 0
						
						neigh_y = int(neigh_id.split('_')[2])
						diag += 1
						
						if neigh_id in _interior_node_ids:
							L[node_ind,neigh_ind] += 1
						elif neigh_id in _exterior_up_node_ids:
							load += self._ly_ind*EY
						elif neigh_id in _exterior_dn_node_ids:
							load += 0
						else:
							raise ValueError('Unknown node id')

						if load != 0:
							B[node_ind,0] -= load
					
						L[node_ind, node_ind] -= diag

				
					# Update conductance 
					C[node_ind, neigh_ind] -= 1

				
				# Record voltage
				self._level_is_solved[_l] = False
				self._level_voltage[_l].append(EY/_max_curr)
				self._level_current[_l].append(_eqns['curr_dn']/_max_curr)
				self._level_scaled_stress[_l].append(_eqns['curr_dn']*ly/(_max_curr*len(self._level_bonds[_l])))
				self._level_stress[_l].append(_eqns['curr_dn']/(lx*_max_curr))
				self._level_strain[_l].append(EY/_max_curr)
				if len(self._level_strain_evp[_l]) == 0:
					self._level_strain_evp[_l].append(EY/_max_curr)
				else:
					self._level_strain_evp[_l].append(EY/_max_curr if EY/_max_curr > self._level_strain_evp[_l][-1] else self._level_strain_evp[_l][-1])
				
			_l -= 1
		self._n_step += 1

	def bond_current(self,bond,l):
		n1_b, n2_b = bond._n1._id, bond._n2._id
		_eqns = self._level_eqns[l]
		curr = _eqns['curr']
		_N_interior, _N_exterior_up, _N_exterior_dn = _eqns['_N_interior'], _eqns['_N_exterior_up'], _eqns['_N_exterior_dn']
		_interior_node_ids, _exterior_up_node_ids, _exterior_dn_node_ids = _eqns['interior_node_ids'], _eqns['exterior_up_node_ids'], _eqns['exterior_dn_node_ids']
		_interior_node_indices, _exterior_up_node_indices, _exterior_dn_node_indices = _eqns['interior_node_indices'], _eqns['exterior_up_node_indices'], _eqns['exterior_dn_node_indices']
		N = _N_interior + _N_exterior_dn + _N_exterior_up
		
		if n1_b in _interior_node_indices:
			n1_ind = _interior_node_indices[n1_b]
		elif n1_b in _exterior_dn_node_indices:
			n1_ind = _exterior_dn_node_indices[n1_b] + _N_interior
		else:
			n1_ind = _exterior_up_node_indices[n1_b] + _N_interior + _N_exterior_dn

		if n2_b in _interior_node_indices:
			n2_ind = _interior_node_indices[n2_b] 
		elif n2_b in _exterior_dn_node_indices:
			n2_ind = _exterior_dn_node_indices[n2_b] + _N_interior
		else:
			n2_ind = _exterior_up_node_indices[n2_b] + _N_interior + _N_exterior_dn

		return curr[n1_ind,n2_ind]


	def simulate_fracture(self,max_step = numpy.inf):
		while True in self._level_not_broken.values() and self._n_step < max_step:
			self.step()
	
	def plot_current(self,**kwargs):
		self._solve()
		_l = self._levels 
		a = self._a

		while _l >= 0:
			_eqns = self._level_eqns[_l]
			curr = _eqns['curr']
			_N_interior, _N_exterior_up, _N_exterior_dn = _eqns['_N_interior'], _eqns['_N_exterior_up'], _eqns['_N_exterior_dn']
			_interior_node_ids, _exterior_up_node_ids, _exterior_dn_node_ids = _eqns['interior_node_ids'], _eqns['exterior_up_node_ids'], _eqns['exterior_dn_node_ids']
			_interior_node_indices, _exterior_up_node_indices, _exterior_dn_node_indices = _eqns['interior_node_indices'], _eqns['exterior_up_node_indices'], _eqns['exterior_dn_node_indices']
			N = _N_interior + _N_exterior_dn + _N_exterior_up
		
			c_mat = numpy.zeros((_eqns['_j_range'],_eqns['_i_range'])) 

			for i in range(N):
				if i < _N_interior:
					node_id = _interior_node_ids[i]
				elif i < _N_interior + _N_exterior_dn:
					node_id = _exterior_dn_node_ids[i-_N_interior]
				else:
					node_id = _exterior_up_node_ids[i-_N_interior-_N_exterior_dn]
				
				for neigh_id in self._level_nodes[_l][node_id]._bonded_to:
					if neigh_id in _interior_node_ids:
						j = _interior_node_indices[neigh_id]
					elif neigh_id in _exterior_dn_node_ids:
						j = _exterior_dn_node_indices[neigh_id] + _N_interior
					elif neigh_id in _exterior_up_node_ids:
						j = _exterior_up_node_indices[neigh_id] + _N_interior + _N_exterior_dn
					else:
						raise ValueError('Unknown neigh_id')

					_ii, _ij = [int(x) for x in node_id.split('_')[1:]]
					_ji, _jj = [int(x) for x in neigh_id.split('_')[1:]]

					#c_mat[(_ii-_eqns['_i_offset']),_ij-_eqns['_j_offset']] += abs(curr[i,j]) 
					#c_mat[(_ji-_eqns['_i_offset']),_jj-_eqns['_j_offset']] += abs(curr[j,i])  
					c_mat[_eqns['_j_range'] - 1 - (_ij-_eqns['_j_offset']),_ii-_eqns['_i_offset']] += abs(curr[i,j]) 
					c_mat[_eqns['_j_range'] - 1 - (_jj-_eqns['_j_offset']),_ji-_eqns['_i_offset']] += abs(curr[j,i])  
				

			fig = mpl.pyplot.figure()
			fig.subplots_adjust(bottom=0.2,left=0.2)
			ax = fig.add_subplot(111)
			ax.set_xscale('linear')
			ax.set_yscale('linear')
			#ax.axis('off')
			if kwargs.get('hide_axis',True):
				ax.set_xticklabels([])
				ax.set_yticklabels([])
				for tic in ax.xaxis.get_major_ticks():
					tic.tick1On = tic.tick2On = False
				for tic in ax.yaxis.get_major_ticks():
					tic.tick1On = tic.tick2On = False
			ax.set_aspect('equal')
	
			#ax.matshow(c_mat.T)
			ax.imshow(c_mat, interpolation='nearest', cmap=cm.jet,alpha=1.0,extent=[0,_eqns['_i_range']*0.5,0,_eqns['_j_range']*(3.0**0.5)*0.5])
			#,extent=[0,cr.cell[0,0],0,cr.cell[1,1]],origin='lower',vmax=kwargs.get('vmax',120),vmin=kwargs.get('vmin',-60))
			mpl.pyplot.draw()

			_l -= 1
	
	def _build_eqns(self):
		"""
		"""

		a = self._a
		ly = self._ly

		_l = self._levels 
		EY = self._EY

		while _l >= 0:
			if not self._level_is_build[_l]:
				# smallest and largest node indices 
				_i_min, _j_min = numpy.inf, numpy.inf
				_i_max, _j_max = -numpy.inf, -numpy.inf

				counter = 0
				_interior_node_ids = []
				_exterior_up_node_ids = []
				_exterior_dn_node_ids = []
				
				for node_id in self._level_nodes[_l]:
					counter += 1
					x, y = [int(z) for z in node_id.split('_')[1:]]
					_i_min = _i_min if _i_min < x else x
					_i_max = _i_max if _i_max > x else x
					_j_min = _j_min if _j_min < y else y
					_j_max = _j_max if _j_max > y else y
				
					if y < self._ly_ind and y > 0:
						_interior_node_ids.append(node_id)
					elif y >= self._ly_ind:
						_exterior_up_node_ids.append(node_id)
					else:
						_exterior_dn_node_ids.append(node_id)
				

				_interior_node_ids.sort()
				_exterior_up_node_ids.sort()
				_exterior_dn_node_ids.sort()

				_N_interior = len(_interior_node_ids)
				_N_exterior_up = len(_exterior_up_node_ids)
				_N_exterior_dn = len(_exterior_dn_node_ids)
				
				_interior_node_indices = {}
				for k in xrange(_N_interior):
					counter+= 1
					node_id = _interior_node_ids[k]
					_interior_node_indices[node_id] = k
			
				_exterior_up_node_indices = {}
				for k in xrange(_N_exterior_up):
					counter+= 1
					node_id = _exterior_up_node_ids[k]
					_exterior_up_node_indices[node_id] = k
				
				_exterior_dn_node_indices = {}
				for k in xrange(_N_exterior_dn):
					counter+= 1
					node_id = _exterior_dn_node_ids[k]
					_exterior_dn_node_indices[node_id] = k
			
				_Li, _Lj, _L_data = [], [], []
				_Bi, _Bj, _B_data = [], [], []
				# Sparse matrix of voltages along the diagonals. Initiated to have interior = -1, exterior up at applied, exterior down at 0
				_Vi, _Vj, _V_data = [], [], []
				# Sparse matrix of conductances
				# Interior, exterior down, exterior up
				_Ci, _Cj, _C_data = [], [], []

				for k in xrange(len(_interior_node_ids)):
					counter+= 1
					node_index = k 
					node_id = _interior_node_ids[k]
					node_x = int(node_id.split('_')[1])
					node_y = int(node_id.split('_')[2])
					neighs = self._level_nodes[_l][node_id]._bonded_to
					diag, load = 0, 0
				
					for neigh_id in neighs:
						counter+= 1
						neigh_y = int(neigh_id.split('_')[2])
						diag += 1
						
						if neigh_id in _interior_node_ids:
							neigh_ind = _interior_node_indices[neigh_id]
							_Li.append(k)
							_Lj.append(neigh_ind)
							_L_data.append(-1.0)

						elif neigh_id in _exterior_up_node_ids:
							load += self._ly_ind*EY

						elif neigh_id in _exterior_dn_node_ids:
							load += 0

						else:
							raise ValueError('Unknown node id')

						if neigh_id in _interior_node_indices:
							neigh_index = _interior_node_indices[neigh_id] 
						elif neigh_id in _exterior_dn_node_indices:
							neigh_index = _exterior_dn_node_indices[neigh_id] + _N_interior
						elif neigh_id in _exterior_up_node_indices:
							neigh_index = _exterior_up_node_indices[neigh_id] + _N_interior + _N_exterior_dn
							
						_Ci.append(node_index)
						_Cj.append(neigh_index)
						_C_data.append(1.0)
					
					if load != 0:
						_Bi.append(k)
						_Bj.append(0)
						_B_data.append(load)
					
					_Li.append(k)
					_Lj.append(k)
					_L_data.append(diag)

					_Vi.append(node_index)
					_Vj.append(node_index)
					_V_data.append(-1.0)
				
				
				
				for k in xrange(len(_exterior_up_node_ids)):
					counter+= 1
					node_index = k + _N_interior + _N_exterior_dn
					node_id = _exterior_up_node_ids[k]
					y = int(node_id.split('_')[2])
					
					_Vi.append(node_index)
					_Vj.append(node_index)
					_V_data.append(y*EY)

				L = scipy.sparse.csc_matrix((_L_data,(_Li,_Lj)),shape=(_N_interior,_N_interior))
				B = scipy.sparse.csc_matrix((_B_data,(_Bi,_Bj)),shape=(_N_interior,1))
				V = scipy.sparse.csc_matrix((_V_data,(_Vi,_Vj)),shape=(_N_interior+_N_exterior_dn+_N_exterior_up,_N_interior+_N_exterior_dn+_N_exterior_up))


				
				for k in xrange(len(_exterior_dn_node_ids)):
					counter+= 1
					node_id = _exterior_dn_node_ids[k]
					node_index = k + _N_interior
					neighs = self._level_nodes[_l][node_id]._bonded_to

					for neigh_id in neighs:
						counter+= 1
						if neigh_id in _interior_node_indices:
							neigh_index = _interior_node_indices[neigh_id] 
						elif neigh_id in _exterior_dn_node_indices:
							neigh_index = _exterior_dn_node_indices[neigh_id] + _N_interior
						elif neigh_id in _exterior_up_node_indices:
							neigh_index = _exterior_up_node_indices[neigh_id] + _N_interior + _N_exterior_dn

						_Ci.append(node_index)
						_Cj.append(neigh_index)
						_C_data.append(1.0)
				
				for k in xrange(len(_exterior_up_node_ids)):
					counter+= 1
					node_id = _exterior_up_node_ids[k]
					node_index = k + _N_interior + _N_exterior_dn
					neighs = self._level_nodes[_l][node_id]._bonded_to

					for neigh_id in neighs:
						counter+= 1
						if neigh_id in _interior_node_indices:
							neigh_index = _interior_node_indices[neigh_id] 
						elif neigh_id in _exterior_dn_node_indices:
							neigh_index = _exterior_dn_node_indices[neigh_id] + _N_interior
						elif neigh_id in _exterior_up_node_indices:
							neigh_index = _exterior_up_node_indices[neigh_id] + _N_interior + _N_exterior_dn

						_Ci.append(node_index)
						_Cj.append(neigh_index)
						_C_data.append(1.0)

				C = scipy.sparse.csc_matrix((_C_data,(_Ci,_Cj)),shape=(_N_interior+_N_exterior_dn+_N_exterior_up,_N_interior+_N_exterior_dn+_N_exterior_up))
				
				ret = {}

				ret['laplacian'] = L
				ret['load'] = B
				ret['conductance'] = C
				ret['conductance_i'] = _Ci
				ret['conductance_j'] = _Cj
				ret['vol_diag'] = V

				ret['interior_node_indices'] = _interior_node_indices
				ret['interior_node_ids'] = _interior_node_ids

				ret['exterior_up_node_ids'] = _exterior_up_node_ids
				ret['exterior_up_node_indices'] = _exterior_up_node_indices
				ret['exterior_dn_node_ids'] = _exterior_dn_node_ids
				ret['exterior_dn_node_indices'] = _exterior_dn_node_indices
				ret['_N_interior'] = _N_interior
				ret['_N_exterior_dn'] = _N_exterior_dn
				ret['_N_exterior_up'] = _N_exterior_up

				ret['_i_range'] = _i_max - _i_min + 1
				ret['_j_range'] = _j_max - _j_min + 1
				ret['_i_offset'] = _i_min
				ret['_j_offset'] = _j_min
				self._level_eqns[_l] = ret
				
				self._level_is_build[_l] = True
			_l -= 1



	def node_coordinates(self,level,boundary_only=False):
		l0 = self._l0
		max_ind_y = self._a*(2*self._ny-1)
		_X, _Y = [], []
		for node_id, node in self._level_nodes[level].iteritems():
			y = int(node_id.split('_')[2])
			if boundary_only:
				if y >= max_ind_y or y <= 0:
					_x, _y = node.coordinates()
					_X.append(_x)
					_Y.append(_y)
			else:
				_x, _y = node.coordinates()
				_X.append(_x)
				_Y.append(_y)

		return _X, _Y


	def plot(self,**kwargs):
		if kwargs.get('show_triangles',False) or kwargs.get('color_by_cur',False):
			self._solve()
		
		cols = ['k','r','b','g']

		"""
		fig = mpl.pyplot.figure()
		fig.subplots_adjust(bottom=0.2,left=0.2)
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')
		#ax.axis('off')
		if kwargs.get('hide_axis',True):
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			for tic in ax.xaxis.get_major_ticks():
				tic.tick1On = tic.tick2On = False
			for tic in ax.yaxis.get_major_ticks():
				tic.tick1On = tic.tick2On = False
		ax.set_aspect('equal')
	
		_l = 0
		while _l <= self._levels:
			if kwargs.get('show_nodes',True):
				x, y = self.node_coordinates(_l)
				line, = ax.plot(x,y,cols[_l]+'.')

			if kwargs.get('show_lines',False):
				for _b_id, _b in self._level_bonds[_l].iteritems():
					_tX,_tY = _b.coordinates(self._nx*self._l0*(self._magnification**(self._levels+1)),self._ny*self._l0*(self._magnification**(self._levels+1)))
					for _X, _Y in zip(_tX,_tY):
						line, = ax.plot(_X,_Y,'b-')

			if kwargs.get('show_notch',False) and self._notch:
				for _b_id, _b in self._level_notched_bonds[_l].iteritems():
					_tX,_tY = _b.coordinates(self._nx*self._l0*(self._magnification**(self._levels+1)),self._ny*self._l0*(self._magnification**(self._levels+1)))
					for _X, _Y in zip(_tX,_tY):
						line, = ax.plot(_X,_Y,'r-')
			
			if kwargs.get('show_broken',True):
				for  _b in self._level_broken_bonds[_l]:
					_tX,_tY = _b.coordinates(self._nx*self._l0*(self._magnification**(self._levels+1)),self._ny*self._l0*(self._magnification**(self._levels+1)))
					for _X, _Y in zip(_tX,_tY):
						line, = ax.plot(_X,_Y,'y-')
			
			_l += 1
		if kwargs.get('show_nodes',True):
			xb, yb = self.node_coordinates(0,boundary_only = True)
			line, = ax.plot(xb,yb,'y.')
		
		mpl.pyplot.draw()
		"""

		_l = 0
		while _l <= self._levels:
			fig = mpl.pyplot.figure()
			#fig.subplots_adjust(bottom=0.2,left=0.2)
			ax = fig.add_subplot(111)
			ax.set_xscale('linear')
			ax.set_yscale('linear')
			#ax.set_title('Level %d'%(_l,))
			ax.axis('off')
			if kwargs.get('hide_axis',True):
				ax.set_xticklabels([])
				ax.set_yticklabels([])
				for tic in ax.xaxis.get_major_ticks():
					tic.tick1On = tic.tick2On = False
				for tic in ax.yaxis.get_major_ticks():
					tic.tick1On = tic.tick2On = False
			ax.set_aspect('equal')

			#ax.set_xlim([0, (2*self._nx - 1)*self._l0*(self._magnification**self._levels)*0.5])		
			#ax.set_ylim([0, (2*self._ny - 1)*self._l0*(self._magnification**self._levels)*(3.0**0.5)*0.5])		
			x, y = self.node_coordinates(_l)
			ax.set_xlim([min(x)-0.1, max(x)+0.1])
			ax.set_ylim([min(y)-0.1, max(y)+0.1])
			if kwargs.get('show_nodes',True):
				x, y = self.node_coordinates(_l)
				line, = ax.plot(x,y,'k.')
			

			if kwargs.get('show_all_triangles',False):
				patches = []
				patch_colors = []
				_temp = self._level_bonds[_l].copy()
				for _b in self._level_broken_bonds[_l]:
					_temp[_b._id] = _b
				_temp.update(self._level_damaged_bonds[_l])


				for _b_id_0, _b_0 in _temp.iteritems():
					i_0, j_0 = [int(x) for x in _b_0._n1._id.split('_')[1:]]
					i_01, j_01 = [int(x) for x in _b_0._n2._id.split('_')[1:]]
					if j_0 == j_01:
						pad = 9
						fmt_str = '_%0' + str(pad) + 'd_%0' + str(pad) + 'd'
						_level_a = self._l0*(self._magnification**_l)
						if i_01 > i_0:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0+2*_level_a,j_0)
							n2_id = fmt_str%(i_0+_level_a,j_0+_level_a)
						else:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0-2*_level_a,j_0)
							n2_id = fmt_str%(i_0-_level_a,j_0+_level_a)


						if False not in [n_id in self._level_nodes[_l] for n_id in [n0_id, n1_id, n2_id]]:
								
							n0 = self._level_nodes[_l][n0_id]
							n1 = self._level_nodes[_l][n1_id]
							n2 = self._level_nodes[_l][n2_id]

							coords = numpy.zeros((3,2))
							for k, n in enumerate([n0_id,n1_id,n2_id]):
								i, j = [int(x) for x in n.split('_')[1:]]
								coords[k,:] = i*0.5*self._l0, j*(3.0**0.5)*0.5*self._l0
						
							patches.append(mpl.patches.Polygon(coords,closed = True))
				
				for _b_id_0, _b_0 in _temp.iteritems():
					i_0, j_0 = [int(x) for x in _b_0._n1._id.split('_')[1:]]
					i_01, j_01 = [int(x) for x in _b_0._n2._id.split('_')[1:]]
					if j_0 == j_01:
						pad = 9
						fmt_str = '_%0' + str(pad) + 'd_%0' + str(pad) + 'd'
						_level_a = self._l0*(self._magnification**_l)
						if i_01 > i_0:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0+2*_level_a,j_0)
							n2_id = fmt_str%(i_0+_level_a,j_0-_level_a)
						else:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0-2*_level_a,j_0)
							n2_id = fmt_str%(i_0-_level_a,j_0-_level_a)
						
						if False not in [n_id in self._level_nodes[_l] for n_id in [n0_id, n1_id, n2_id]]:
						
							n0 = self._level_nodes[_l][n0_id]
							n1 = self._level_nodes[_l][n1_id]
							n2 = self._level_nodes[_l][n2_id]

							coords = numpy.zeros((3,2))
							for k, n in enumerate([n0_id,n1_id,n2_id]):
								i, j = [int(x) for x in n.split('_')[1:]]
								coords[k,:] = i*0.5*self._l0, j*(3.0**0.5)*0.5*self._l0
						
							patches.append(mpl.patches.Polygon(coords,closed = True))
				p = mpl.collections.PatchCollection(patches,facecolors='None',edgecolors='blue',alpha=0.5, lw=0.5)
				ax.add_collection(p)
			
			
			if kwargs.get('show_triangles',False):
				patches = []
				patch_colors = []

				c_min = abs(self._level_eqns[_l]['curr'].data).min()
				c_max = abs(self._level_eqns[_l]['curr'].data).max()

				map_cur = lambda x : (x - c_min)/(c_max - c_min)

				for _b_id_0, _b_0 in self._level_bonds[_l].iteritems():
					i_0, j_0 = [int(x) for x in _b_0._n1._id.split('_')[1:]]
					i_01, j_01 = [int(x) for x in _b_0._n2._id.split('_')[1:]]
					if j_0 == j_01:
						pad = 9
						fmt_str = '_%0' + str(pad) + 'd_%0' + str(pad) + 'd'
						_level_a = self._l0*(self._magnification**_l)
						if i_01 > i_0:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0+2*_level_a,j_0)
							n2_id = fmt_str%(i_0+_level_a,j_0+_level_a)
						else:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0-2*_level_a,j_0)
							n2_id = fmt_str%(i_0-_level_a,j_0+_level_a)


						if False not in [n_id in self._level_nodes[_l] for n_id in [n0_id, n1_id, n2_id]]:
								
							n0 = self._level_nodes[_l][n0_id]
							n1 = self._level_nodes[_l][n1_id]
							n2 = self._level_nodes[_l][n2_id]


							if n0.is_bonded_to(n1) and n0.is_bonded_to(n2) and n1.is_bonded_to(n2):
								if n0_id + '_' + n1_id in self._level_bonds[_l]:
									b01 = self._level_bonds[_l][n0_id + '_' + n1_id]
								else:
									b01 = self._level_bonds[_l][n1_id + '_' + n0_id]
								
								if n0_id + '_' + n2_id in self._level_bonds[_l]:
									b02 = self._level_bonds[_l][n0_id + '_' + n2_id]
								else:
									b02 = self._level_bonds[_l][n2_id + '_' + n0_id]
								
								if n1_id + '_' + n2_id in self._level_bonds[_l]:
									b12 = self._level_bonds[_l][n1_id + '_' + n2_id]
								else:
									b12 = self._level_bonds[_l][n2_id + '_' + n1_id]

								cur = max( [abs(self.bond_current(b,_l)) for b in [b01, b02, b12]])
								coords = numpy.zeros((3,2))
								for k, n in enumerate([n0_id,n1_id,n2_id]):
									i, j = [int(x) for x in n.split('_')[1:]]
									coords[k,:] = i*0.5*self._l0, j*(3.0**0.5)*0.5*self._l0
							
								color = cm.jet(map_cur(cur))
								patches.append(mpl.patches.Polygon(coords,closed = True, facecolor=color, edgecolor='none'))
								patch_colors.append(color)
				
				for _b_id_0, _b_0 in self._level_bonds[_l].iteritems():
					i_0, j_0 = [int(x) for x in _b_0._n1._id.split('_')[1:]]
					i_01, j_01 = [int(x) for x in _b_0._n2._id.split('_')[1:]]
					if j_0 == j_01:
						pad = 9
						fmt_str = '_%0' + str(pad) + 'd_%0' + str(pad) + 'd'
						_level_a = self._l0*(self._magnification**_l)
						if i_01 > i_0:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0+2*_level_a,j_0)
							n2_id = fmt_str%(i_0+_level_a,j_0-_level_a)
						else:
							n0_id = fmt_str%(i_0  ,j_0)
							n1_id = fmt_str%(i_0-2*_level_a,j_0)
							n2_id = fmt_str%(i_0-_level_a,j_0-_level_a)
						
						if False not in [n_id in self._level_nodes[_l] for n_id in [n0_id, n1_id, n2_id]]:
						
							n0 = self._level_nodes[_l][n0_id]
							n1 = self._level_nodes[_l][n1_id]
							n2 = self._level_nodes[_l][n2_id]


							if n0.is_bonded_to(n1) and n0.is_bonded_to(n2) and n1.is_bonded_to(n2):
								if n0_id + '_' + n1_id in self._level_bonds[_l]:
									b01 = self._level_bonds[_l][n0_id + '_' + n1_id]
								else:
									b01 = self._level_bonds[_l][n1_id + '_' + n0_id]
								
								if n0_id + '_' + n2_id in self._level_bonds[_l]:
									b02 = self._level_bonds[_l][n0_id + '_' + n2_id]
								else:
									b02 = self._level_bonds[_l][n2_id + '_' + n0_id]
								
								if n1_id + '_' + n2_id in self._level_bonds[_l]:
									b12 = self._level_bonds[_l][n1_id + '_' + n2_id]
								else:
									b12 = self._level_bonds[_l][n2_id + '_' + n1_id]

								cur = sum( [abs(self.bond_current(b,_l)) for b in [b01, b02, b12]])/3
								cur = max( [abs(self.bond_current(b,_l)) for b in [b01, b02, b12]])
								coords = numpy.zeros((3,2))
								for k, n in enumerate([n0_id,n1_id,n2_id]):
									i, j = [int(x) for x in n.split('_')[1:]]
									coords[k,:] = i*0.5*self._l0, j*(3.0**0.5)*0.5*self._l0
							
								color = cm.jet(map_cur(cur))
								patches.append(mpl.patches.Polygon(coords,closed = True, facecolor=color, edgecolor='none'))
								patch_colors.append(color)
				p = mpl.collections.PatchCollection(patches,facecolors=patch_colors,edgecolors=patch_colors)
				ax.add_collection(p)

			if kwargs.get('show_notch',False) and self._notch:
				for _b_id, _b in self._level_notched_bonds[_l].iteritems():
					_tX,_tY = _b.coordinates(self._nx*self._l0*(self._magnification**(self._levels+1)),self._ny*self._l0*(self._magnification**(self._levels+1)))
					for _X, _Y in zip(_tX,_tY):
						line, = ax.plot(_X,_Y,'r-')
			
			if kwargs.get('show_damage',False) and self._damage:
				for _b_id, _b in self._level_damaged_bonds[_l].iteritems():
					_tX,_tY = _b.coordinates(self._nx*self._l0*(self._magnification**(self._levels+1)),self._ny*self._l0*(self._magnification**(self._levels+1)))
					for _X, _Y in zip(_tX,_tY):
						line, = ax.plot(_X,_Y,'k-',lw=2)
			
			if kwargs.get('show_broken',True):
				for  _b in self._level_broken_bonds[_l]:
					_tX,_tY = _b.coordinates(self._nx*self._l0*(self._magnification**(self._levels+1)),self._ny*self._l0*(self._magnification**(self._levels+1)))
					for _X, _Y in zip(_tX,_tY):
						line, = ax.plot(_X,_Y,'r-',lw=2)
			
			if kwargs.get('show_nodes',True):
				xb, yb = self.node_coordinates(_l,boundary_only = True)
				line, = ax.plot(xb,yb,'y.')
			
			if kwargs.get('show_lines',False):
				if kwargs.get('color_by_cur',False):
					c_min = abs(self._level_eqns[_l]['curr'].data).min()
					c_max = abs(self._level_eqns[_l]['curr'].data).max()

					map_cur = lambda x : (x - c_min)/(c_max - c_min)

				for _b_id, _b in self._level_bonds[_l].iteritems():
					_tX,_tY = _b.coordinates(self._nx*self._l0*(self._magnification**(self._levels+1)),self._ny*self._l0*(self._magnification**(self._levels+1)))
					for _X, _Y in zip(_tX,_tY):
						_X, _Y = numpy.array(_X), numpy.array(_Y)
						if kwargs.get('color_by_cur',False):
							color = cm.jet(map_cur(abs(self.bond_current(_b,_l))))
							line, = ax.plot(_X,_Y,color=color)
						else:
							line, = ax.plot(_X,_Y,'b-',alpha=0.5,lw=0.5)
							
			mpl.pyplot.draw()
			if self._notch:
				#fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f.eps'%(self._nx, self._ny,_l,self._magnification,self._level_notch_len[_l]),format='eps',pad_inches=0,bbox_inches='tight',dpi=144)
				fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f.png'%(self._nx, self._ny,_l,self._magnification,self._level_notch_len[_l]),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
				fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f.pdf'%(self._nx, self._ny,_l,self._magnification,self._level_notch_len[_l]),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
			else:
				#fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f.eps'%(self._nx, self._ny,_l,self._magnification,0),format='eps',pad_inches=0,bbox_inches='tight',dpi=144)
				fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f.png'%(self._nx, self._ny,_l,self._magnification,0),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
				fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f.pdf'%(self._nx, self._ny,_l,self._magnification,0),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
			
			_l += 1


def scan_matrix(m):
	# return the max element and its i,j index in the compressed sparse column matrix m 
	_col = 0
	for k in xrange(len(m.data)):
		if k >= m.indptr[_col+1]:
			_col += 1
		if abs(m.data[k] - m[m.indices[k],_col]) < 1E-3:
			print [k, m.indices[k], _col]
		else:
			raise Exception('%d, %d, %d, %.2f, %.2f'%(k, m.indices[k], _col, m.data[k], m[m.indices[k],_col]))

def csc_argmax(m, I=None, J=None):
	if I is None:
		# return the max element and its i,j index in the compressed sparse column matrix m 
		_max = -numpy.inf 
		_col = 0
		_max_data_ind = 0
		_max_col = 0
		for k in xrange(len(m.data)):
			if k >= m.indptr[_col+1]:
				_col += 1
			if m.data[k] > _max:
				_max = m.data[k]
				_max_data_ind = k
				_max_col = _col
		
		return _max, m.indices[_max_data_ind], _max_col
	
	else:
		_max = -numpy.inf 
		_i_max, _j_max = 0, 0
		for i, j in zip(I,J):
			if m[i,j] > _max:
				_i_max, _j_max = i, j
				_max = m[i,j]

		return _max, i, j
			
		

def csc_argmin(m):
	# return the max element and its i,j index in the compressed sparse column matrix m 
	_min = numpy.inf 
	_col = 0
	_min_data_ind = 0
	_min_col = 0
	for k in xrange(len(m.data)):
		if k >= m.indptr[_col+1]:
			_col += 1
		if m.data[k] < _min:
			_min = m.data[k]
			_min_data_ind = k
			_min_col = _col

	return _min, m.indices[_min_data_ind], _min_col


def sample():
	levels = [0,1]
	level_size = {0:[64], 1:[10]}
	level_damage = {0: [0.01,0.05,0.1,0.15,0.2,0.3], 1: [0.01,0.05,0.1,0.15,0.2,0.3]}
	level_sims = {0:100, 1:100}

	res = {'levels': levels, 'level_size': level_size, 'level_damage': level_damage, 'level_sims': level_sims}
	for lev in levels:
		lev_res = {}
		for lev_size in level_size[lev]:
			for damage in level_damage[lev]:
				avg_max_stress = 0
				for sim in range(level_sims[lev]):
					print sim
					simulated = False
					while not simulated:
						try:
							hg = hierarchical_grid(lev_size,lev_size,lev,l0=1,notch=False,damage=True,damage_fraction=damage)
							hg.simulate_fracture()
							simulated = True
						except ValueError:
							simulated = False
					lev_res['stress'] = hg._level_stress[0]
					lev_res['strain'] = hg._level_strain[0]
					avg_max_stress += max(lev_res['stress'])
					res['level_%d_size_%d_damage_%.2f_sim_%d'%(lev,lev_size,damage,sim)] = lev_res
				avg_max_stress /= level_sims[lev]
				res['avg_max_stress_level_%d_size_%d_damage_%.2f'%(lev,lev_size,damage)] = avg_max_stress

	return res

def stress(hg,l):
	_i_min, _j_min = numpy.inf, numpy.inf
	_i_max, _j_max = -numpy.inf, -numpy.inf
	
	for node_id in hg._level_nodes[l]:
		x, y = [int(z) for z in node_id.split('_')[1:]]
		_i_min = _i_min if _i_min < x else x
		_i_max = _i_max if _i_max > x else x
		_j_min = _j_min if _j_min < y else y
		_j_max = _j_max if _j_max > y else y

	# FIXME
	_notch_len_ind = 0.25*(_i_max - _i_min)
	
	_j_mid = 0.5*(_j_min + _j_max)
	if abs( round(_j_mid) - _j_mid)  < 1E-3:
		_j_mid += 0.25


	dist = []
	curr = []


	for _b_id, _b in hg._level_bonds[l].iteritems():
		n1_id, n2_id = _b._n1._id, _b._n2._id
		n1_id, n2_id = _b._n1._id, _b._n2._id
		n1_x, n1_y = [int(x) for x in n1_id.split('_')[1:]]
		n2_x, n2_y = [int(x) for x in n2_id.split('_')[1:]]
		if n1_x - _i_min >= _notch_len_ind and n2_x - _i_min >= _notch_len_ind and (n1_y - _j_mid)*(n2_y - _j_mid) < 0:
			dist.append( ( n1_x + n2_x - 2*_notch_len_ind)/2.0)
			curr.append(abs(hg.bond_current(_b,l)))

	dist = numpy.array(dist)
	curr = numpy.array(curr)
	sort_ind = numpy.argsort(dist)
	dist = dist[sort_ind]
	curr = curr[sort_ind]

	return dist*0.5, curr

def t_stress(c,d):
	return c[0] + c[1]*(d**-0.5)

def t_stress_diff(c,d,s):
	return s - t_stress(c,d)

def fit_stress(d, s):
	c = scipy.optimize.leastsq(t_stress_diff,numpy.zeros(2),(d,s))[0]
	return c, t_stress(c,d)


def rt_inv_diff(c,x,y):
	return y - c[0]*(x**-0.5)

def fit_rt_inv(x, y):
	c = scipy.optimize.leastsq(rt_inv_diff,numpy.zeros(1),(x,y))[0]
	return c, c*(x**-0.5)

def rt_diff(c,x,y,w):
	return y - c[0]*(x**0.5) - w*c[1]*(x**-0.5)

def fit_rt(x, y, w):
	c = scipy.optimize.leastsq(rt_diff,numpy.zeros(2),(x,y,w))[0]
	return c, c[0]*(x**0.5) + w*c[1]*(x**-0.5)

def den_diff(c,x,y):
	return y - c[0]/x - c[1]/(x**2.0)

def fit_den(x, y):
	c = scipy.optimize.leastsq(den_diff,numpy.zeros(2),(x,y))[0]
	return c, c[0]/x + c[1]/(x**2.0)

#-----
def fit_str_dist_func(c,m,d,stress, surv):
	y = numpy.log(-numpy.log(surv)/(17*10))
	x = c[0] + (c[1] + c[3]/m) * (5.0*c[2]*c[2]*numpy.log(d) + numpy.log(m))/(stress*stress*m*m)
	return x, y
	
def fit_str_dist_diff(c,data):
	err = numpy.array([])
	for k, v in data.iteritems():
		d, m = [float(x) for x in k.split('_')]
		stress, surv = v
		x, y = fit_str_dist_func(c,m,d,stress,surv)
		err = numpy.concatenate((err, y - x))
	
	return err

def fit_str_dist(data):
	c = scipy.optimize.leastsq(fit_str_dist_diff,numpy.array([0.12, 0.05, 1.0, 3.0]),(data,))[0]
	return c


#-----
def fit_str_func(c,mag,damage, string=False, ft='raw'):
	if ft == 'raw':
		if string:
			return '%.2f/m * ( %.2f log(p) - log(m) )**0.5'%(c[3],-5*(c[2]**2.0))
		else:
			return (c[3]/mag) * (-5*(1.0*c[2]**2.0)*numpy.log(damage) - numpy.log(mag))**0.5

	elif ft == 'den':
		if string:
			return '(%.2f + %.2f/m) * ( %.2f log(p) - log(m) )**0.5'%(c[0], c[1], -5*(c[2]**2.0))
		else:
			return (c[0] + c[1] / (mag**1.0)) * (-5*(1.0*c[2]**2.0)*numpy.log(damage) - numpy.log(mag))**0.5 
	
def str_diff(c,mag,damage,strength, strength_den, std_strength, std_strength_den):
	err = numpy.array([])
	for d in damage:
		y = strength[d]
		wt = std_strength[d]
		err = numpy.concatenate((err, (y - fit_str_func(c,mag, d, ft='raw'))/(y) ))
		y = strength_den[d]
		wt = std_strength_den[d]
		err = numpy.concatenate((err, (y - fit_str_func(c,mag, d, ft='den'))/(y) ))
	return err 

def fit_str(mag, damage, strength, strength_den, std_strength, std_strength_den):
	c = scipy.optimize.leastsq(str_diff,numpy.array([0.12, 0.05, 2.0, 1.0]),(mag, damage, strength, strength_den, std_strength, std_strength_den))[0]
	return c

#----------

def plot_stress(hg,l=0,x0=-1):
	if x0 == -1:
		x0 = hg._magnification/8

	x, s = stress(hg,l)
	x += x0 - x[0]
	max_cur = max(abs(hg._level_eqns[l]['curr'].data))
	s = s/max_cur
	s = s/hg._level_density[l]

	s_avg = numpy.convolve(s, numpy.ones((5,))/5, mode='valid')
	x = x[:len(x)-4]
	s = s[:len(s)-4]

	par, f = fit_stress(x, s_avg)

	fig = mpl.pyplot.figure()
	ax = fig.add_subplot(111)
	ax.set_xscale('linear')
	ax.set_yscale('linear')
	ax.set_xlim([0, max(x)])

	ax.plot(x,s_avg,'ko',alpha=0.5)
	ax.plot(x,f,'k-')
	ax.plot(x,s,'ro-',markeredgecolor='r',alpha=0.5)
	[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
	[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
	fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.png'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
	fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.pdf'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)

	
	fig = mpl.pyplot.figure()
	ax = fig.add_subplot(111)
	ax.set_xscale('linear')
	ax.set_yscale('linear')
	ax.set_xlim([60, 80])
	ax.set_ylim([0.1, 0.25])

	ax.plot(x,s_avg,'ko',alpha=0.5)
	ax.plot(x,f,'k-')
	ax.plot(x,s,'ro-',markeredgecolor='r',alpha=0.5)
	[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
	[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
	fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress_Zoom.png'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
	fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress_Zoom.pdf'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)

	
def toughness_stress(mag=None, toughness_fit=None,toughness_notch=None,density=None, levels = 1):

	"""
	if mag is None:
		mag = [8, 16, 24, 32, 40, 48, 56, 64]
		toughness_fit = []
		toughness_notch = []
		density = []
		for m in mag:
			print m
			hg0 = hierarchical_grid(17,10,levels,magnification=m,l0=1,notch=True)
			hg0._solve()
			x, s = stress(hg0,0)
			x += m/8 - x[0]
			max_cur = max(abs(hg0._level_eqns[0]['curr'].data))
			s = s/max_cur
			s = s*hg0._level_density[0]
			s_avg = numpy.convolve(s, numpy.ones((5,))/5, mode='valid')
			x = x[:len(x)-4]
			s = s[:len(s)-4]

			par, f = fit_stress(x, s_avg)
			toughness_fit.append(par[1])
			toughness_notch.append( (hg0._level_eqns[0]['curr_dn'] / (max_cur*hg0._lx))*(hg0._level_notch_len[0]**0.5))
			density.append(hg0._level_density[0])

		res = {}
		res['mag'] = mag
		res['toughness_fit'] = toughness_fit
		res['toughness_notch'] = toughness_notch
		res['density'] = density
		res['x'] = x
		res['stress'] = s
		f_out = open('temp.pic','w+')
		pickle.dump(res,f_out)
		f_out.close()
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')
		ax.set_title('Magnification = %d'%(m))
		ax.set_xlim([0, max(x)])

		ax.plot(x,s_avg,'ko',alpha=0.5)
		ax.plot(x,f,'k-')
		ax.plot(x,s,'ro-',markeredgecolor='r',alpha=0.5)
		[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
		mpl.pyplot.draw()
		fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.png'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
		fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.pdf'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)

	with open('Nx_17_Ny_10_NotchData.pic','r') as f:
		data = pickle.load(f)

	mag = numpy.array(mag)
	toughness_fit	= numpy.array(toughness_fit)
	toughness_notch	= numpy.array(toughness_notch)
	density = numpy.array(density)
	
	"""
	with open('temp.pic','r') as f:
		data = pickle.load(f)
	
	mag = data['mag']
	toughness_fit = data['toughness_fit']
	density = data['density']
	
	mag = numpy.array(mag)
	toughness_fit	= numpy.array(toughness_fit)
	density = numpy.array(density)

	fig = mpl.pyplot.figure()
	ax = fig.add_subplot(111)
	ax.set_xscale('linear')
	ax.set_yscale('linear')
	#ax.set_xlim([0, max(x)])

	print mag
	print toughness_fit 
	c, f = fit_rt_inv(mag, toughness_fit)

	ax.plot(mag,toughness_fit,'ko')
	x = numpy.linspace(mag.min(), mag.max(),100)
	ax.plot(x,c[0]*x**-0.5,'k-')
	print "K = %.2f /sqrt(m)"%(c[0])

	#ax.set_yticks([0.04,0.08,0.12, 0.16])
	ax.set_xticks([0,8,16,32,64])
	[tick.label.set_fontsize(30) for tick in ax.yaxis.get_major_ticks()]
	[tick.label.set_fontsize(30) for tick in ax.xaxis.get_major_ticks()]
	mpl.pyplot.draw()
	fig.savefig('Fit_K_Mag_Level_%d.png'%(levels),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
	fig.savefig('Fit_K_Mag_Level_%d.pdf'%(levels),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
	
	
	fig = mpl.pyplot.figure()
	ax = fig.add_subplot(111)
	ax.set_xscale('linear')
	ax.set_yscale('linear')
	#ax.set_xlim([0, max(x)])
	
	tough_per_density = toughness_fit/density
	c, f = fit_rt(mag, tough_per_density, 1.0)

	ax.plot(mag,tough_per_density,'ko')
	x = numpy.linspace(mag.min(), mag.max(),100)
	ax.plot(x,c[0]*x**0.5 + c[1]*x**-0.5,'k-')
	print "K/rho = %.2f x**0.5 + %.2f/x**0.5"%(c[0],c[1])
	
	c, f = fit_rt(mag, tough_per_density, 0.0)
	#ax.plot(x,c[0]*x**0.5 + 0*c[1]*x**-0.5,'k--')
	print "K/rho = %.2f x**0.5 + 0*%.2f/x**0.5"%(c[0],c[1])
	ax.set_yticks([0.4,0.6,0.8, 1.0])
	ax.set_xticks([0,8,16,32,64])
	[tick.label.set_fontsize(30) for tick in ax.yaxis.get_major_ticks()]
	[tick.label.set_fontsize(30) for tick in ax.xaxis.get_major_ticks()]
	
	mpl.pyplot.draw()
	fig.savefig('Fit_K_Den_Mag_Level_%d.png'%(levels),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
	fig.savefig('Fit_K_Den_Mag_Level_%d.pdf'%(levels),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)

	"""
	fig = mpl.pyplot.figure()
	ax = fig.add_subplot(111)
	ax.set_xscale('linear')
	ax.set_yscale('linear')
	#ax.set_xlim([0, max(x)])
	
	c, f = fit_rt_inv(mag, toughness_notch)

	ax.plot(mag,toughness_notch,'ko',alpha=0.5)
	x = numpy.linspace(mag.min(), mag.max(),100)
	ax.plot(x,c[0]*x**-0.5,'k-')
	[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
	[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
	mpl.pyplot.draw()
	fig.savefig('Notch_K_Mag_Levels_%d.png'%(levels),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
	fig.savefig('Notch_K_Mag_Levels_%d.pdf'%(levels),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
	"""
	

	return mag, toughness_fit, toughness_notch, density	


def strength_plots():
	
	mag = numpy.array([8, 16, 32, 64])
	density = numpy.array([2.422603626548774, 1.3686602635865615,0.7236697443713243, 0.37166977533017304])
	damage_f = numpy.array([0.1, 0.2])

	cols = ['k','r','b']


	if 1 == 2:
		hg = hierarchical_grid(17,10,1,magnification=8,l0=1,notch=False,damage=True,damage_fraction=0.1)
		hg.simulate_fracture()
		
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')
		ax.plot(hg._level_stress[0],'ko-')
		[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
		fig.savefig('StressCurve.png',format='png',pad_inches=0,bbox_inches='tight',dpi=144)
		fig.savefig('StressCurve.pdf',format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
		
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')
		ax.plot(hg._level_strain[0],'ko-')
		[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
		fig.savefig('StrainCurve.png',format='png',pad_inches=0,bbox_inches='tight',dpi=144)
		fig.savefig('StrainCurve.pdf',format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
	
	
	if 1 == 1:
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')

		# Gather data
		all_data = {}
		for i, df in enumerate(damage_f):
			for m, den in zip(mag,density):
				f = open('Nx_17_Ny_10_Mag_%d_StrengthData.pic'%(m),'r')
				data = pickle.load(f)
				f.close()
				stress = data['mag_%d_damage_%.2f_stress'%(m,df)]
				stress.sort()
				surv = 1.0 - (0.35 + numpy.arange(len(stress)))/len(stress)
				all_data['%.2f_%d'%(df, m)] = (stress, surv)
		
		# Fit 
		c = fit_str_dist(all_data)
		print c
		for i, df in enumerate(damage_f):
			for m, den in zip(mag,density):
				stress, surv = all_data['%.2f_%d'%(df,m)]
				x, y = fit_str_dist_func(c,m,df,stress,surv)
				ax.plot(x,y,'k.')
		#ax.set_yticks([0.04,0.08,0.12, 0.16])
		#ax.set_xticks([0,8,16,32,64])
		[tick.label.set_fontsize(30) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(30) for tick in ax.xaxis.get_major_ticks()]
	
	if 1 == 2:
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')

		# Gather data
		mean_strength = {}
		mean_strength_den = {}
		std_strength = {}
		std_strength_den = {}
		for i, df in enumerate(damage_f):
			mean_str_den = []
			mean_str = []
			std_str_den = []
			std_str = []
			for m, den in zip(mag,density):
				f = open('Nx_17_Ny_10_Mag_%d_StrengthData.pic'%(m),'r')
				data = pickle.load(f)
				f.close()
				mean_str_den.append(data['mag_%d_damage_%.2f_stress'%(m,df)].mean()/den)
				mean_str.append(data['mag_%d_damage_%.2f_stress'%(m,df)].mean())
				
				std_str_den.append(data['mag_%d_damage_%.2f_stress'%(m,df)].std()/den)
				std_str.append(data['mag_%d_damage_%.2f_stress'%(m,df)].std())
			
			mean_str_den = numpy.array(mean_str_den)
			mean_strength_den[df] = mean_str_den
			std_str_den = numpy.array(std_str_den)
			std_strength_den[df] = std_str_den
			
			mean_str = numpy.array(mean_str)
			mean_strength[df] = mean_str
			std_str = numpy.array(std_str)
			std_strength[df] = std_str

		# Fit 
		c = fit_str(mag, damage_f, mean_strength, mean_strength_den, std_strength, std_strength_den)
		x = numpy.linspace(mag[0], mag[-1], 100)
		for i, df in enumerate(damage_f):
			mean_str_den = mean_strength_den[df]
			std_str_den = std_strength_den[df]
			print mean_str_den
			print std_str_den
			ax.errorbar(mag,mean_str_den,yerr=std_str_den*1.5, fmt=cols[i]+'o')
			f = fit_str_func(c, x, df, ft='den') 	
			ax.plot(x,f,cols[i]+'-')
		print fit_str_func(c,0,0,string=True, ft='den')
		ax.set_yticks([0.04,0.08,0.12, 0.16])
		ax.set_xticks([0,8,16,32,64])
		[tick.label.set_fontsize(30) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(30) for tick in ax.xaxis.get_major_ticks()]
		
		fig.savefig('Strength_PerDen.png',format='png',pad_inches=0,bbox_inches='tight',dpi=144)
		fig.savefig('Strength_PerDen.pdf',format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
	
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')

		for i, df in enumerate(damage_f):
			mean_str = mean_strength[df]
			std_str = std_strength[df]
			ax.errorbar(mag,mean_str,yerr=std_str*1.5,fmt=cols[i]+'o')
			f = fit_str_func(c, x, df, ft='raw') 	
			ax.plot(x,f,cols[i]+'-')
		print fit_str_func(c,0,0,string=True, ft='raw')
		ax.set_yticks([0,0.1,0.2,0.3,0.4])
		ax.set_xticks([0,8,16,32,64])
		[tick.label.set_fontsize(30) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(30) for tick in ax.xaxis.get_major_ticks()]
		fig.savefig('Strength_Raw.png',format='png',pad_inches=0,bbox_inches='tight',dpi=144)
		fig.savefig('Strength_Raw.pdf',format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
	
	
	if 1 ==2:
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')

		ax.plot(mag,density,'ko')
		c, f = fit_den(mag,density)
		x = numpy.linspace(mag.min(), mag.max(), 100)
		f = c[0]/x + c[1]/(x*x)
		ax.plot(x,f,'k-')
				
		[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
		#fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.png'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
		#fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.pdf'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
	
	if 1 == 2:
		fig = mpl.pyplot.figure()
		ax = fig.add_subplot(111)
		ax.set_xscale('linear')
		ax.set_yscale('linear')

		for i, df in enumerate(damage_f):
			mean_str = []
			for m, den in zip(mag,density):
				print (m, df)
				f = open('Nx_17_Ny_10_Mag_%d_StrengthData.pic'%(m),'r')
				data = pickle.load(f)
				f.close()
				stress = data['mag_%d_damage_%.2f_stress'%(m,df)]
				stress.sort()
				_cdf = (numpy.arange(len(stress)) + 0.35)/len(stress)
				_sdf = 1.0 - _cdf
				ax.plot((1.0 / (m * m)) * ((1.0 / stress)**2.0) * (-7 * numpy.log(df) - numpy.log(m)), numpy.log(-numpy.log(_sdf)/(m*m)), 'k.')
				#ax.plot((1.0/m)*(1.0/stress)**2.0, numpy.log(-numpy.log(_sdf)), 'k.')
		[tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
		[tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
		#fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.png'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='png',pad_inches=0,bbox_inches='tight',dpi=144)
		#fig.savefig('Nx_%d_Ny_%d_Lev_%d_Mag_%d_Notch_%.2f_Stress.pdf'%(hg._nx, hg._ny,l,hg._magnification,hg._level_notch_len[l]),format='pdf',pad_inches=0,bbox_inches='tight',dpi=144)
	
		
