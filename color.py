'''
Library of miscellaneous functions related to color handling.
'''

from numpy import cos, pi, sqrt, log10
from matplotlib.colors import LinearSegmentedColormap
from func_lib import *
from colorist import ColorHex
from typing import Literal
import matplotlib.pyplot as plt

def procHex(*args):
	'''
	Converts the given values (0-255) into a hex color code.
	'''
	return '#'+''.join([f'0{hex(s)[2:]}'[-2:] for s in args])

def str2Color(s: str):
	'''
	Converts the given string into an arbitrary (but consistent) color code.
	'''
	tot_num = 3*2551*sum([ord(c)**5 for c in s])
	R = int(tot_num%255)
	G = int((tot_num/1739)%255)
	B = int((tot_num*3717)%255)
	if R + G + B > 2.5*255: [R, G, B] = [int(c/2) for c in [R, G, B]]
	return procHex(R, G, B)

def pop2Alpha(p: str):
	'''
	Returns an alpha value corresponding with the kind of population (S, I, R) submitted. (R is lighter.)

	*Largely obsolete.*
	'''
	p = p[0]
	if p == 'S' or p == 'I': return 1.
	else: return 0.7

# Basic functions for color mapping. Determines how each channel's value changes over a 0-1 interval
def r(i): return sqrt(int(-127*cos(pi*i) + 127)/255)
def g(i): return int(-200*cos(0.5*pi*i) + 200)/255
def b(i): return int(70*cos(0.75*pi*(0.5*i + 0.25)) + 70)/255

def getColor(rv: float, gv: float, bv: float):
	'''
	Turns the given RGB values (floats, range 0-1) into a hexadecimal color code.
	'''
	return procHex(*[int(rgb*255) for rgb in [r(rv), g(gv), b(bv)]])

def scale_mid(lsp: list[float]):
	'''
	Scales the given list (should be linear, bounded between 0 and 1) such that values near the middle (0.5) are "flattened." Useful
	for color maps intending to exaggerate differences near the edges (0, 1).
	'''
	return [((2*i-1)**3 + 1)/2 for i in lsp]

def scale_edg(lsp: list[float]) -> list[float]:
	'''
	Scales the given list (should be linear, bounded between 0 and 1) such that values near the edges (0, 1) are "flattened." Useful
	for color maps intending to exaggerate differences near the middle (0.5).
	'''
	return [(np.sign(2*i-1)*log10(1 + 9*abs(2*i-1)) + 1)/2 for i in lsp]

def scale_lin(lsp: list[float]):
	'''
	"Scales" the given list (should be linear, bounded between 0 and 1) such that its values are linearly separated, i.e. leaves them
	unchanged. Mainly a function wrapper indicating that no action should be taken.
	'''
	return lsp

def getColorMap(nm: str='cmap', res: int=256, scale_type: Literal['mid','edge','lin'] = 'lin'):
	'''
	Creates a `mpl.LinearSegmentedColormap` object based on the gradient defined by `r`, `g`, and `b` (dark blue-orange-gold). `res` is
	its resolution, `nm` is its name, and `scale_type` defines the method used to scale the color map's values.
	- `mid`: Flattens the middle (0.5), exaggerates the edges (0, 1).
	- `edge`: Flattens the edges (0, 1), exaggerates the middle (0.5).
	- `lin`: Leaves the scaling unchanged (default).
	'''
	scale_fn = scale_mid
	if scale_type is 'edge': scale_fn = scale_edg
	elif scale_type is 'lin': scale_fn = scale_lin
	elif scale_type is not 'mid': print('invalid scale type submitted; defaulting to mid')
	def makeGrad(fn):
		idxs = linspace(res=res)
		idxs_unsc = scale_fn(idxs)
		r_t = ()
		for i in range(res):
			idx = idxs[i]
			idx_u = idxs_unsc[i]
			y1 = fn(idxs_unsc[i-1])
			if not i: y1 = fn(idx_u)
			r_t += ((idx, fn(idx_u), y1),)
		return r_t

	return LinearSegmentedColormap(nm, {'red': makeGrad(r), 'green': makeGrad(g), 'blue': makeGrad(b)})

def plotGradient():
	'''
	Displays the gradient defined by `r`, `g`, and `b` (dark blue-orange-gold).
	'''
	idxs = linspace()
	test_spc = linspace(res=5)
	hexes = [getColor(i, i, i) for i in test_spc]
	print(', '.join(f'{ColorHex(hex_col)}{hex_col[1:]}{ColorHex(hex_col).OFF}' for hex_col in hexes))
	plt.imshow((idxs,idxs), aspect='auto', cmap=getColorMap())
	plt.show()