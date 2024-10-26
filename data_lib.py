'''
Library of miscellaneous functions related to data analysis and recording.
'''

from func_lib import *
from typing import TextIO
import pandas as pd

def propUnc(*args, do_div: bool=True) -> float:
	'''
	Propagates uncertainty for mean calculations. Takes an arbitrary number of standard deviation arguments. `do_div` denotes whether or not
	this is for a mean's direct uncertainty (`True`) or a sum of a standard deviation with an uncertainty (`False`).
	'''
	div_fac = 1
	if do_div: div_fac = len(args)
	return np.sqrt(sum([(std/div_fac)**2 for std in args]))

def writeOptLine(f: TextIO, name: str, mean_num: float, std_num: float):
	'''
	Writes the given data (`mean_num`, `std_num`) to file `f` with row name `name`.
	'''
	f.write(f'{name}:\t{roundAndSN(mean_num)}\t+- {roundAndSN(std_num)}\n')

def readCSVData(csv_path: str) -> dict[str, list[float]]:
	'''
	Returns data from a CSV file as a dict of column header:value list.
	'''
	return {col_header: list(col_data) for col_header, col_data in pd.read_csv(csv_path).items()}

def saveStats(lst: list[float], frac_to_take: float=0.2) -> tuple[float, float]:
	'''
	Gets the mean & standard deviation of the last `frac_to_take` proportion of its elements. Returned in that order.
	'''
	data_to_keep = lst[int(-frac_to_take*len(lst)):]
	return np.mean(data_to_keep), np.std(data_to_keep)

def getRange(base_val: float, dim_val: float, num_els: int=5):
	'''
	Gets a range around `base_val` with dimensions `dim_val` in each direction & total element count `num_els`.
	'''
	u_lim = base_val + dim_val
	l_lim = base_val - dim_val
	step_size = (u_lim - l_lim)/(num_els-1)
	return [l_lim + i*step_size for i in range(num_els)]

def makeCoordDF(c_dict: dict[tuple[float, float], float]):
	'''
	Turns the given dictionary (structure - coordinate tuple:value, i.e. (x, y):z) into a `pandas.DataFrame` object.
	'''
	df_data = [(y, x, z) for (x, y), z in c_dict.items()]
	return pd.DataFrame(df_data, columns=['y', 'x', 'z']).pivot(index='y', columns='x', values='z')