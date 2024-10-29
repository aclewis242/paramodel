'''
Library of miscellaneous functions, most of which being fairly general and not semantically tied to a specific project.
'''

from typing import Any
import numpy as np
import os

def normalise(l: list[float]):
	'''
	Normalises the given list such that its sum is 1.
	'''
	l_sum = sum(l)
	if not l_sum: return l
	return [l_el/l_sum for l_el in l]

def normalise_np(l: np.ndarray[float]):
	'''
	Normalises the given NumPy array.
	'''
	return l/sum(l)

def normalise_dict(d: dict[str, float]) -> dict[str, float]:
	'''
	Normalises the given dictionary (as per its values).
	'''
	return dictify(d.keys(), normalise(d.values()))

def fn(f_n: str):
	'''
	Converts the given string into a proper file name.
	'''
	return f_n.lower().replace(' ', '_')

def float2SN(f: float, p: int=2, do_sci: bool=False):
	'''
	Converts the given number to scientific notation with the given precision.

	### Parameters
	- `f`: The number in question.
	- `p`: The number (int) of digits of precision (past the decimal point). For example: f=1302, p=2 -> 1.30e3.
	- `do_sci`: Whether or not to force scientific notation (e.g. 13 -> 1.3e1).
	'''
	if ((f < 10**(p+1) and f > 10**(-p)) or not f) and not do_sci: return f'{f}'
	else:
		pwr = np.floor(np.log10(f))
		return f'{int(f/(10**(pwr-p)))/(10**p)}e{int(pwr)}'

def roundNum(f: float, prec: int=2) -> float:
	'''
	Rounds the given number to the given number of decimal points.

	### Parameters
	- `f`: The number in question.
	- `prec`: The number of decimal points to round it to.
	'''
	return round(f*(10**prec))/(10**prec)

def roundAndSN(f: float, u_lim: int=4, l_lim: int=2, prec: int=3):
	'''
	General-purpose number processing method. Will return a number with the given precision if it falls outside the given limits.

	### Parameters
	- `f`: The number in question.
	- `u_lim`: The power of 10 to use as an upper limit for scientific notation.
	- `l_lim`: The power of 10 to use as a lower limit for scientific notation. Positive means negative (i.e. l_lim = 2 means a limit of 0.01.)
	'''
	if (f < 10**u_lim and f > 10**(-l_lim+1)) or not f: return f'{roundNum(f, prec=prec)}'
	else: return float2SN(f, p=prec, do_sci=True)

def roundTrailing(*args, max_prec: int=6) -> list[float]:
	'''
	Rounds trailing 0s/9s that may arise from tiny floating-point errors elsewhere.
	'''
	rv = []
	for arg in args:
		arg_str = str(arg)
		if len(arg_str) > max_prec + 2:
			if not int(arg_str[-2]): rv += [float(arg_str[:-1])]
			else: rv += [arg + 10**(2-len(arg_str))]
		else: rv += [arg]
	return rv

def normPercentList(l: list[float]) -> list[float]:
	'''
	Turns the given list into a normalised list of percentages.
	'''
	return list(100*normalise_np(np.array(l)))

def printFloatList(l: list[float]):
	'''
	Prints the items of a list of floats (rounded to 2 decimal places), as well as their indices.
	'''
	[print(f'{i}:\t{roundNum(l[i])}') for i in range(len(l))]

def printMat(m: list[list]):
	'''
	Prints the given 2d matrix (usually of transition probabilities, though not necessarily) in an easier-to-read fashion.
	'''
	rv = '\n'.join(['\t'.join([str(int(100*i)) for i in j]) for j in m])
	print(rv)
	return rv

def hapify(g: str):
	'''
	Turns the given genotype into its haploid equivalent (e.g., AA.Bb.cc becomes A.B.c). Haploid genotypes are unaffected.

	*Note that heterozygous genes are not preserved, and are instead assumed effectively mutated/capital.*
	*As heterozygotes comprise a minuscule fraction of the overall population, the effect of this is negligible.*
	'''
	return '.'.join([s[0] for s in g.split('.')])

def dipify(g: str):
	'''
	Turns the given genotype into its diploid equivalent (e.g., A.B.c becomes AA.BB.cc). Diploid genotypes are unaffected.
	'''
	return '.'.join([''.join(2*[s])[:2] for s in g.split('.')])

def listify(a: list | np.ndarray | Any) -> list | Any:
	'''
	Produces a deep copy of type list from a given multi-dimensional list (or NumPy array). If the elements of the list are complex
	types, then these objects will still be the same in memory (the deep copy only extends to the lists themselves).
	'''
	if (type(a) is np.ndarray) or (type(a) is list): return [listify(el) for el in a]
	else: return a

def dictify(ks: list, vs: list):
	'''
	Turns the given lists of keys and values into a dict. These should be the same length!
	'''
	return {k: v for k, v in zip(ks, vs)}

def mkDir(*args):
	'''
	Makes a new directory with the given name, if it does not exist already. Takes an arbitrary number of arguments.
	'''
	[os.mkdir(arg) for arg in args if arg not in os.listdir()]

def mkFilePath(*args):
	'''
	Makes a new file at the given path, if it does not exist already. Takes an arbitrary number of arguments.
	'''
	for a in args:
		if a not in os.listdir(): f = open(a, 'x'); f.close()

def mkFile(f_path: str):
	'''
	Makes and returns a new file at the given path, clearing it if it already exists.
	'''
	if os.path.exists(f_path): os.remove(f_path)
	return open(f_path, 'x')

def getSubDirs(f_path: str) -> list[str]:
	'''
	Returns a list of all directories within the given directory (does not include full path or their own sub-directories).
	'''
	if f_path[-1] != '/': f_path += '/'
	if not os.path.exists(f_path): return []
	return [path_loc for path_loc in os.listdir(f_path) if os.path.isdir(f'{f_path}{path_loc}')]

def numSubDirs(f_path: str):
	return len(getSubDirs(f_path)) if os.path.exists(f_path) else 0

def transpose(l: list[list]):
	'''
	Transposes the given two-dimensional list. Must be rectangular; that is, all the second-order lists must be equal in length.
	'''
	return [list(l2) for l2 in zip(*l)]

def linspace(l_lim: float=0., u_lim: float=1., res: int=256):
	'''
	Generates an evenly-spaced list of `res` points between `l_lim` and `u_lim`.
	'''
	return [l_lim + u_lim*i/(res-1) for i in range(res)]

def isInRange(num: float, rng: list[float]=[0.,1.], do_trunc_check: bool=True):
	'''
	Checks whether or not the given number is in the given range (default 0-1).
	'''
	if do_trunc_check: num = float(trunc(num)); rng = [float(trunc(r)) for r in rng]
	return num >= rng[0] and num <= rng[1]

def trunc(pop_val: str | float, trunc_len: int=7, **kwargs):
	'''
	Truncates the given string/float (cast to string) such that it has the specified length (default 7).
	'''
	pop_val = str(pop_val)
	return pop_val if len(pop_val) <= trunc_len else pop_val[:trunc_len]

def list_str(lst: list, limit: int=40, shoulder: int=5):
	'''
	Writes the given list as a string, if it is too long to reasonably print to the console.
	'''
	if len(lst) <= limit: return str(lst)
	else: return f'{str(lst[:shoulder])[:-1]} ... {str(lst[-shoulder:])[1:]}'