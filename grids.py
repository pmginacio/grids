# -*- coding: utf-8; -*-
'''A library to handle gridded data in python.'''
import numpy as np
import utils
import cPickle as pickle

__author__ = "Pedro Inácio"
__copyright__ = "Copyright 2020"
__version__ = "1.0"
__maintainer__ = "Pedro Inácio"
__email__ = "pedromiragaia@gmail.com"
__status__ = "Development"

class Grid(object):
	"""class definition for data equi-spaced in a 2D grid"""
	def __init__(self, *args, **kwargs):
		"""
		A grid may be created:
		A) empty from an x, y array of coordinates
			A1) Grid(x=x,y=y) TODO
			A2) Grid(xlim=[xmin, xmax],ylim=[ymin,ymax],dx=dx,dy=dy)
		B) from an existing array of values, with or without coordinates
		C) from a complete list or partial list of x,y,z triplets, each row
			with one triplet
			C1) Grid(xyz=xyz, dx=None, dy=None, xlim=None, ylim=None)
		"""
		super(Grid, self).__init__()

		# get input x and y lists
		xyz = kwargs.pop('xyz',None)
		if xyz:
			xin = xyz[:,0]
			yin = xyz[:,1]
			zin = xyz[:,2]

			# compute lims and spacing
			xin.sort(), yin.sort()
			xlimin = [xin[0], xin[-1]]
			ylimin = [yin[0], yin[-1]]

			# find minimum spacing between coordinate values
			dxin = np.min(xin[1:] - xin[0:-1])
			dyin = np.min(yin[1:] - yin[0:-1])

		# get optional grid parameters and check consistency with xyz input
		self.xlim = kwargs.pop('xlim', xlimin)
		if self.xlim[0] > xlimin[0] or self.xlim[1] < xlimin[1]:
			raise ValueError('xlim bounds are smaller than the input xyz list of points')

		self.ylim = kwargs.pop('ylim', ylimin)
		if self.ylim[0] > ylimin[0] or self.ylim[1] < ylimin[1]:
			raise ValueError('ylim bounds are smaller than the input xyz list of points')

		self.dx = kwargs.pop('dx',dxin)
		if self.dx > dxin:
			raise ValueError('dx is larger than minimin spacing in xyz list of points')

		self.dy = kwargs.pop('dy',dyin)
		if self.dy > dyin:
			raise ValueError('dy is larger than minimin spacing in xyz list of points')

		# generate grid coordinates
		# compute number of grid nodes
		self.nx = (self.xlim[1] - self.xlim[0]) / float(self.dx) + 1
		self.ny = (self.ylim[1] - self.ylim[0]) / float(self.dy) + 1
		# nx should be very close to an integer
		if not utils.num_equal(self.nx, int(self.nx)):
			raise ValueError('dx does not divide the xyz list in an integer number of nodes')
		if not utils.num_equal(self.ny, int(self.ny)):
			raise ValueError('dy does not divide the xyz list in an integer number of nodes')

		self.x = np.linspace(self.xlim[0],self.xlim[1],self.nx)
		self.y = np.linspace(self.ylim[0],self.ylim[1],self.ny)

		# check that all xin and yin are present in self.x and self.y, otherwise,
		# the input list of points is not regularly spaced.
		if not all([utils.num_equal(np.min(self.x - i), 0.0) for i in xin]):
			raise ValueError('x list of points is not regularly spaced')
		if not all([utils.num_equal(np.min(self.y - i), 0.0) for i in yin]):
			raise ValueError('y list of points is not regularly spaced')

		# arrays are stored with y-along the rows and x- along the columns
		self.z = np.array((self.ny,self.nx),np.nan)
		# store input list
		for x,y,z in xyz:
			ix = int(x-self.xlim[0]/float(dx))
			iy = int(y-self.ylim[0]/float(dy))
			self.z[i,j] = z

		# optional metadata
		self.xlabel = kwargs.pop('xlabel', None)
		self.xunit = kwargs.pop('xunit', None)

		self.ylabel = kwargs.pop('ylabel', None)
		self.yunit = kwargs.pop('yunit', None)

		self.zlabel = kwargs.pop('zlabel', None)
		self.zunit = kwargs.pop('zunit', None)
	

	# def pack(self):

	# def unpack(z):

	def save(self, file):
		'''
		save the grid to a file using pickle
		'''

		with open(file, 'wb') as f:
			#  pickle using the highest protocol available.
			pickle.dump(self, f, -1)

	def regrid(self, grd):

	def plot(self, file):
		import matplotlib.pyplot as plt

		fig, ax = plt.subplots()
		im = ax.imshow(trend, origin='lower', cmap='bwr', vmax=10, vmin=-10)
		fig.colorbar(im, ax=ax, orientation='horizontal', fraction=.1)
		fig.savefig(file)

	def __str__(self):
		'''
		Print information about the grid
		'''

		aux = 'Grid:\n'
		aux += 'Shape:	%dx%d\n' % (self.nx, self.ny)
		aux += 'Limits:	%f-%f,%f-%f\n' % (self.xlim, self.ylim)
		aux += 'Spacing:	%f,%f\n' % (self.dx, self.dy)
		aux += 'Units:	%s, %s, %s\n' % (self.xunit, self.yunit, self.zunit)
		aux += 'Label:	%s, %s, %s\n' % (self.xlabel, self.ylabel, self.zlabel)
		aux += 'Data:\n%s\n' % (str(self.data))

	# implement operators
	# article that covers the subject
	#	https://realpython.com/operator-function-overloading/

def read(file):
	'''load grid from the pickle file'''

	with open(file, 'rb') as f:
		grd = pickle.load(file)

	return grd