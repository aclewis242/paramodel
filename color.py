from numpy import cos, pi
from matplotlib.colors import LinearSegmentedColormap

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
	'''
	p = p[0]
	if p == 'S' or p == 'I': return 1.
	else: return 0.7

def getColorMap(nm: str='cmap', res: int=256):
	'''
	Calculates RGB color based on position in a custom gradient and returns it as a hexadecimal color code.

	### Parameters
	i: A decimal between 0 and 1 representing position within the gradient.
	'''
	def r(i): return int(-120*cos(pi*i) + 120)/256
	def g(i): return int(-120*cos(0.5*pi*i) + 120)/256
	def b(i): return int(70*cos(0.75*pi*(0.5*i + 0.25)) + 70)/256
	def makeGrad(fn):
		idxs = [i/(res-1) for i in range(res)]
		r_t = ()
		for i in range(res):
			idx = idxs[i]
			y1 = fn(idxs[i-1])
			if not i: y1 = fn(idx)
			r_t += ((idx, fn(idx), y1),)
		return r_t

	return LinearSegmentedColormap(nm, {'red': makeGrad(r), 'green': makeGrad(g), 'blue': makeGrad(b)})