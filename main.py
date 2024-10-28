'''
Main file for the host-vector parasite & SIRS model.
'''

from sim_lib import *
from allele import *
from color import *
from data_lib import *
from time import time
from scipy.interpolate import griddata
from scipy.stats import linregress
from PIL import Image, ImageDraw
from io import BytesIO
import matplotlib.pyplot as plt
import pandas as pd
import os
import cProfile
import pstats
import colorist

H1 = 'h1'
VEC = 'vec'
FULL_NMS = {'Vector': VEC, 'Host': H1}
FULL_NMS_R = {VEC: 'Vector', H1: 'Host'}

# Rates are in terms of '# events expected per day'
VCT = {				# Parameters for transmission vector behavior (mosquito)
	'bd': 0.071,	# Birth/death rate
	'ir': 0.,		# Infection rate between individuals. Largely obsolete for a vector-based model
	'rr': 0,		# Recovery rate
	'wi': 0.,		# Waning immunity rate
	'pn': VEC,				# Short population name, used in console output & within the code
	'pn_full': 'Vector',	# Full population name, used in graph titles
	'is_vector': True,		# Whether or not the population is a disease vector
}

hst_transm_p_wi = 0.11		# The base transmission probability for hosts, only used here to come up with wi_rate
ct_rate = 0.28				# Contact rate (average bites per day)
tau = 90					# Average duration of immunity (days)
h_ir = ct_rate*hst_transm_p_wi								# Infection rate (only used in the context of wi_rate)
wi_rate = h_ir*np.exp(-h_ir*tau)/(1 - np.exp(-h_ir*tau))	# Waning immunity rate

HST1 = {
	'bd': 0.,
	'ir': 0.,
	'rr': 2.19e-3,
	'wi': wi_rate,
	'pn': H1,
	'pn_full': 'Host',
}

para_lifespan = 8.				# Parasite lifespan (hours)
INDV_VEC = {
	'pc': 20,					# Parasite count
	'para_lsp': para_lifespan, 
	'is_hap': False,			# Whether or not this individual's parasites are haploid
	'do_mutation': False,		# Whether or not this individual's parasites can mutate
	'pc_to_transmit': 15,		# The number of parasites transmitted per infection
}
INDV_HST = {
	'pc': int(1e8),
	'mut_chance': 2.94e-6,		# Chance of mutation per parasite per generation
	'para_lsp': para_lifespan,
	'is_hap': True,
	'do_mutation': True,
	'pc_to_transmit': 10,
}
INDVS = [INDV_VEC, INDV_HST]

D = allele(char='D')

####################################
### ----- BEGIN USER INPUT ----- ###
####################################

# Selection biases for mutated allele, relative to the wild allele (always 1.0)
control_sa = {H1: 1.0, VEC: 1.0}		# Used for the control case
hst_sa_inv = {H1: 1.05, VEC: 1.0}		# Used for the invasion example of the host selection advantage case
hst_sa_exp = {H1: 1.002, VEC: 1.0}		# Used for the host selection advantage case (quantitative)
vec_sa_exp = {H1: 1.0, VEC: 1.5}		# Used for the vector selection advantage case
asel_explc = {H1: 1.002, VEC: 0.5}		# Used for the antagonistic selection explicit case

asel_dims = {H1: 0.004, VEC: 0.5}		# Range to use in each direction around the base value for antagonistic

D.sel_advs = control_sa		# Choose the parameter set to use here (except antagonistic, which is further below)

# For transmission probabilities: the pop ID is the source -- i.e., vec: 0.021 means 2.1% transmission chance from vector to host

hst_base_transm_p = 0.11	# Probability of contact resulting in a host-vector transmission
hst_adv_transm_p = 0.12		# Advantageous transmission probability for hosts

vec_base_transm_p = 0.021	# Probability of contact resulting in a vector-host transmission
vec_adv_transm_p = 0.024	# Advantageous transmission probability for vectors

D.base_transm_probs = {H1: hst_base_transm_p, VEC: vec_base_transm_p} # for wild-type allele

base_ta = D.base_transm_probs.copy()					# No advantages
hst_ta = {H1: hst_adv_transm_p, VEC: vec_base_transm_p} # Host transmission advantage
vec_ta = {H1: hst_base_transm_p, VEC: vec_adv_transm_p} # Vector transmission advantage
atrans = {H1: 0.10, VEC: vec_adv_transm_p}				# Antagonistic transmission

atrans_dims = {H1: 0.01, VEC: 0.003}	# Range in each direction to use around the base value

D.transm_probs = atrans 	# Choose the parameter set to use here (except antagonistic, which is below)

antag_types = {'sel': asel_dims, 'trans': atrans_dims}

NUM_RUNS = 1				# Number of simulations to run
FILE_DIR = 'hta' 	# Directory to save files under, if multiple simulations are being run (irrelevant for antagonistic)
INIT_MUT_PROP = 0.5			# Initial proportion of mutated alleles
SIM_LENGTH = 20000			# The length of the simulation (days)
SHOW_RES = True 			# Whether or not to show the results
PROP_FOR_STATS = 0.2		# The proportion (from the end) of the results to use for statistics
NUM_HISTS = 0				# Set to 0 for no histograms

DO_ANTAG = True 			# Whether or not to do antagonistic parameters
ANTAG_TYPE = 'trans'  		# 'sel' | 'trans' (for antagonistic selection and transmission respectively)
CONTOUR_DENSITY = 5 		# Performs (num)^2 rounds
DO_RETRO_CONTOUR = False 	# Whether or not to generate the contour plot from (all) existing data. Will not generate anything new
DO_EXTEND = True			# Whether or not to extend preexisting data. Will either generate or read data, depending on what's available
COLOR_SCALE_TYPE = 'edge'	# 'lin' | 'mid' | 'edge', how to scale the colors. 'mid' flattens the middle and 'edge' flattens the edges

##################################
### ----- END USER INPUT ----- ###
##################################

TO_RUN = 'primary'	# 'cost' | 'primary' | 'overall'

ANTAG_NMS = {'sel': 'selection', 'trans': 'transmission'}

USER_INPTS = {'Simulation time': SIM_LENGTH,
			  'Total number of runs': NUM_RUNS,
			  'Selection advantages': D.sel_advs,
			  'Transmission advantages': D.transm_probs}

SHOW_CONTOUR = SHOW_RES
SHOW_MULTI = SHOW_RES and not DO_ANTAG
if DO_ANTAG or NUM_RUNS > 1: SHOW_RES = False

ALLELES = [D] # Do NOT have more than one allele here -- the simulation has been optimised for the single-locus case.
			  # Adding more WILL break it!

PARAMS_1 = HST1
PARAMS_2 = VCT

def run(p0: np.ndarray=np.array([[20, 1, 0], [21, 0, 0]], dtype='float64'), p_fac: float=1200., nt: float=1., num_hist: int=NUM_HISTS,
		plot_res: bool=SHOW_RES, t_scale: float=SIM_LENGTH, init_mut_prop: float=INIT_MUT_PROP, fdir: str=''):
	'''
	Run the simulation. Keyword arguments should be changed in the "user input" section of the main file (above this method).

	### Parameters
	- `p0`: The initial population ratios (S, I, R) as a NumPy array of 3-element NumPy arrays.
	- `p_fac`: The scale factor on population.
	- `nt`: Time steps per day. Alternatively, 24/`nt` is hours per time step.
	- `num_hist`: The number of histograms to generate (varying over time).
	- `plot_res`: Whether or not to display the results' graph, as a bool.
	- `t_scale`: The number of days to run the simulation for.
	- `init_mut_prop`: The initial fraction of mutated (uppercase) alleles. Must be between 0 and 1.
	- `fdir`: The directory to save the results under.

	### Returns
	- `data_to_return`: dict of infected population name to a two-element list containing the mean (0) and standard deviation (1) of the data.
	'''
	exts_to_rm = ['dat', 'csv', 'txt']
	[[os.remove(file) for file in os.listdir() if file.endswith(f'.{ext}')] for ext in exts_to_rm]
	mkDir('hists', 'old images')
	if fdir: os.makedirs(fdir, exist_ok=True)
	mkFilePath('inf_events_raw.dat', 'last_event.dat',)
	alleles = ALLELES
	p0 = p_fac*p0
	t_max = t_scale
	for i in range(len(INDVS)):
		i_params = INDVS[i]
		para_gens = int((24/nt)/i_params['para_lsp'])
		i_params['para_gens'] = para_gens
		i_params['alleles'] = alleles
	nt = float(int(nt*t_scale))

	hosts_1 = population(p0[0], **INDV_HST)
	vectors = population(p0[1], **INDV_VEC)

	m1 = SIR(hosts_1, **PARAMS_1)
	m2 = SIR(vectors, **PARAMS_2)

	m1.itr = {vectors: ct_rate}
	m1.other_pop = vectors
	m2.itr = {hosts_1: ct_rate}
	m2.other_pop = hosts_1

	mdls = [m1, m2]

	t0 = time()
	ts, ps, pops, hists_v, hists_h, hist_tms, mdls_fin = simShell(
		t_max, mdls, nt=nt, alleles=alleles, init_mut_prop=init_mut_prop, num_hist=num_hist)
	ex_tm = time() - t0
	print(f'\nExecution time: {roundNum(ex_tm, prec=3)}') # consider colored console output for readability

	[p.printDat() for p in pops]
	[mdl.save(f'{fdir}models') for mdl in mdls_fin]
	# dimensions of ps: layer 1 is times, layer 2 is models at that time, layer 3 is pop #s for that model

	f = open('inf_events.dat', 'x')
	f_raw = open('inf_events_raw.dat', 'r')
	inf_events = {}
	for line in f_raw.readlines():
		line = line[:-1]
		if line in inf_events: inf_events[line] += 1
		else: inf_events[line] = 1
	[f.write(f'{i_e}: {inf_events[i_e]}\n') for i_e in inf_events]
	f.close()
	f_raw.close()

	def getDims(lst: list, tab: str=''): # Displays the dimensions of the given list. Useful when handling complex/unorthodox structures
		if type(lst) == list:
			print(f'{tab}list of dim {len(lst)} containing:')
			getDims(lst[0], f'{tab}\t')
			if type(lst[0]) == list:
				if len(lst[1]) != len(lst[0]):
					print(f'{tab}additionally:')
					getDims(lst[1], f'{tab}\t')
			return
		else:
			print(f'{tab}{type(lst).__name__}')
			return
	
	output_fn = f'{fdir}net_output.opt' # '.opt' is a plain text file. Only marked that way to keep the file clearing from catching it
	f = mkFile(output_fn)
	data_to_return: dict[str, list[float]] = {}
	for i in range(len(mdls)):
		ns = [''.join(n.split('.')) for n in pops[i].getAllPopNms()]
		f.write(f'\t{mdls[i].pn_full}:\n')
		gens = []
		for n in ns:
			if '(' in n and ')' in n: gens += [n[n.index('(')+1:n.index(')')]]
			else: gens += [n]
		ps_i = np.array([k[i] for k in ps])
		csv_data = {'times': ts}
		plt_datas = []
		stplt_labels = []
		stplt_colors = []
		for j in range(len(ns)):
			if (ns[j][0] != 'R') and (ns[j][0] != 'S'):
				plt_data = ps_i[:,j]
				plt_lst = list(plt_data)
				plt_datas += [plt_lst]
				stplt_labels += [ns[j]]
				plt_color = str2Color(gens[j])
				stplt_colors += [plt_color]
				csv_data[ns[j]] = plt_lst
				mean_stdev = saveStats(plt_lst, frac_to_take=PROP_FOR_STATS)
				data_to_return[ns[j]] = list(mean_stdev)
				writeOptLine(f, ns[j], *mean_stdev)
				plt.plot(ts, plt_data, label=ns[j], color=plt_color, alpha=pop2Alpha(ns[j]))
		f.write('\n')

		# Plot frequencies
		net_i = [1. for t in ts]
		plt.plot(ts, net_i, label='I (total)')
		plt.plot(ts, 0*ts, alpha=0.)
		plt.title(f'{mdls[i].pn_full} population (infected)')
		plt.legend()
		plt.xlabel('Simulation time')
		plt.ylabel('Population')
		file_nm = f'{fdir}{fn(mdls[i].pn)}'
		plt.close()

		pd.DataFrame(csv_data).to_csv(f'{file_nm}.csv')

		# Plot proportions
		plt.stackplot(ts, *plt_datas, labels=stplt_labels, colors=stplt_colors)
		plt.title(f'{mdls[i].pn_full} population: proportion plot')
		plt.legend()
		plt.ylabel('Relative proportion')
		plt.xlabel('Simulation time')
		plt.savefig(f'{file_nm}.png')
		if plot_res: plt.show()
		plt.close()
	
	# Save histograms
	hists_with_pn = {VEC: hists_v, H1: hists_h}
	hist_max = p0[0][1]
	for hpn in hists_with_pn:
		hists_2_p = hists_with_pn[hpn]
		for i in range(num_hist):
			tm = hist_tms[i]
			hist_2_p = [0.0] + hists_2_p[i] + [1.0]
			plt.hist(hist_2_p, bins=100)
			plt.title(f'Population {hpn}, time {roundNum(tm)}')
			plt.xlabel('Mutated (capital) allele frequencies')
			plt.ylabel('Individual count')
			if hpn == H1: plt.ylim(top=hist_max)
			plt.savefig(f'hists/{hpn}_{i}.png')
			plt.close()
	
	writeInputs(f)
	f.close()
	return data_to_return

def doTimeBreakdown(): # Runs the simulation with a profiler & processes the output accordingly
	time_output_fn = 'time_breakdown.dat'
	cProfile.run('run()', time_output_fn)
	time_output_txt = 'time_breakdown.txt'
	time_output_file = open(time_output_txt, 'x')
	p_stats = pstats.Stats(time_output_fn, stream=time_output_file)
	p_stats.strip_dirs()
	p_stats.sort_stats('cumtime')
	p_stats.print_stats()

def doOverallData(all_data: dict[str, list[list[float]]], full_dir: str=''):
	'''
	Takes the results of multiple simulations, writes them to a file, and returns them.

	### Parameters
	- `all_data`: dict of strain (I (D), I (Dd), I (d), ...):list of mean-standard deviation lists (at indices 0 and 1 respectively).
	   Each element should be the result of a different simulation.
	- `full_dir`: The directory to save the output file under (full filepath).

	### Returns
	- `f`: The file object the overall data was saved to.
	- `stat_lsts`: The data saved in `f` (dict strain:[mean, standard deviation]).
	'''
	all_data: dict[str, list[list[float]]] = {k: transpose(all_data[k]) for k in all_data}
	lengths = {'Vector': 6, 'Host': 5} # lengths of row names (ref. net_output.opt et al for examples)
	f = mkFile(f'{full_dir}net_output_overall.opt')
	f.write(f'\t----- OVERALL OUTPUT: {full_dir.split("/")[-2]} -----\n')
	stat_lsts: dict[str, list[float]] = {k: [np.mean(v[0]), propUnc(np.std(v[0]), propUnc(*v[1]), do_div=False)] for k, v in all_data.items()}
	for pn, lng in lengths.items():
		keys = [k for k in all_data if len(k) == lng]
		keys.sort(reverse=True)
		f.write(f'\t{pn}:\n')
		[writeOptLine(f, k, *stat_lsts[k]) for k in keys]
		f.write('\n')
	writeInputs(f)
	return f, stat_lsts

def writeInputs(f: TextIO):
	[f.write(f'{k}: {v}\n') for k, v in USER_INPTS.items()]

def doMultipleRuns(n: int=3, fdir: str='', force_multi: bool=False, do_qc: bool=False) -> dict[str, list[float]]:
	'''
	Performs multiple simulations, saves the data, and returns the results (dict strain:[mean, standard deviation]).

	### Parameters
	- `n`: The number of simulations to perform.
	- `fdir`: The directory (under full_outputs) to save the output files to.
	- `force_multi`: Whether or not to force the multiple-run structure to be used. If `False` and `n` is 1, it will simply perform a single
	   run as normal.
	'''
	if n == 1 and not force_multi: run(); return {}
	full_dir = f'full_outputs/{fdir}/'
	all_data: dict[str, list[list[float]]] = {}
	num_runs_done = numSubDirs(full_dir)
	num_tot = n + num_runs_done
	ignore_amt = 0.1
	for i in range(1,num_tot+1):
		if num_tot != 1: colorist.red(10*'-' + f' RUN {i} ' + 10*'-')
		ret_data = {}
		if i <= num_runs_done:
			colorist.red(f'\tRun {i} already complete. Skipping.')
			ret_data = doRetroStats(fdir)
		else: ret_data = run(fdir=f'{full_dir}{i}/')
		for k in ret_data:
			if k in all_data: all_data[k] += [ret_data[k]]
			else: all_data[k] = [ret_data[k]]
		if do_qc:
			final_freq = getNetAllFreq(ret_data)
			if not isInRange(final_freq, [ignore_amt, 1-ignore_amt], do_trunc_check=False):
				colorist.red(f'\tResult insufficiently ambiguous to warrant further simulations. Breaking.')
				break
	f, stat_lsts = doOverallData(all_data, full_dir)
	f.close()
	compilePlots(full_dir[:-1], do_show=SHOW_MULTI)
	return stat_lsts

def doRetroStats(fdir: str):
	'''
	Retroactively calculates overall stats from the given directory (under full_outputs). Useful if data has been generated, but the
	`stat_lsts` objects from those generations are not available for whatever reason (this method rebuilds them).
	'''
	full_dir = f'full_outputs/{fdir}/'
	csv_names = FULL_NMS
	bad_cols = ['Unnamed', 'times']
	all_data: dict[str, list[list[float]]] = {}
	for r_dir in os.listdir(full_dir):
		if '.' in r_dir: continue
		full_r_dir = f'{full_dir}{r_dir}/'
		if not os.listdir(full_r_dir): continue
		for csv_nm in csv_names.values():
			csv_data: dict[str, list[float]] = readCSVData(f'{full_r_dir}{csv_nm}.csv')
			real_bad_cols: list[str] = []
			for k in csv_data:
				for bad_col in bad_cols:
					if bad_col in k: real_bad_cols += [k]
			for rbc in real_bad_cols: del csv_data[rbc]
			for k, v in csv_data.items():
				stats_to_add = list(saveStats(v, frac_to_take=PROP_FOR_STATS))
				if k in all_data: all_data[k] += [stats_to_add]
				else: all_data[k] = [stats_to_add]
	f, stat_lsts = doOverallData(all_data, full_dir)
	f.close()
	return stat_lsts

def getRetroContourData(fn_save: str):
	'''
	Reads data from the given file name and returns it in the format required of the contour plotting methods.

	### Parameters
	- `fn_save`: The full path to the desired file, not including the .csv extension.

	### Returns
	- `xs_vec`: The column names of the csv file. Corresponds to vector parameters, which are on the x-axis.
	- `ys_h1`: The row names of the csv file. Corresponds to host parameters, which are on the y-axis.
	- `zs_freq`: A two-dimensional list of the values in the csv file, structured similarly.
	'''
	contour_df = pd.read_csv(f'{fn_save}.csv')
	xs_vec = [float(x) for x in contour_df.columns if x != 'y']
	ys_h1 = [float(y) for y in contour_df.loc[:,'y']]
	zs_freq: list[list[float]] = [[ys[i_x] for i_x in range(len(ys)) if i_x] for ys in contour_df.to_numpy()]
	return xs_vec, ys_h1, zs_freq

def getGenContourData(pop_data: dict[tuple[float,float], float], fn_save: str):
	'''
	Saves the contour plot data to a particular file and returns it in the format required of the plotting methods.

	### Parameters
	- `pop_data`: A dict of coordinates to frequencies (x,y):z used to generate a contour plot.
	- `fn_save`: The full path to the desired file, not including the .csv extension.

	### Returns
	- `xs_vec`: The column names of the csv file. Corresponds to vector parameters, which are on the x-axis.
	- `ys_h1`: The row names of the csv file. Corresponds to host parameters, which are on the y-axis.
	- `zs_freq`: A two-dimensional NumPy array of the values in the csv file, structured similarly.
	'''
	contour_df = makeCoordDF(pop_data)
	contour_df.to_csv(f'{fn_save}.csv')
	xs_vec = [float(x) for x in contour_df.columns]
	ys_h1 = [float(y) for y in contour_df.index]
	zs_freq = contour_df.to_numpy()
	return xs_vec, ys_h1, zs_freq

def getNetAllFreq(stat_lsts: dict[str, list[float]], pop_len: Literal[5,6]=5):
	'''
	Takes a `stat_lsts` object (dict strain:[mean, standard deviation]), filters out the relevant entries (`pop_len` 5 for hosts, 6 for
	vectors), and returns the net frequency of the mutated allele across all strains. Does not consider standard deviation.
	'''
	keys = [k for k in stat_lsts if len(k) == pop_len]
	keys.sort(reverse=True)
	net_all_freq = 0.
	for k in keys: net_all_freq += k.count(D.char)*stat_lsts[k][0]
	net_all_freq /= (len(keys) - 1)
	return net_all_freq

def doContourPlots(retro: bool=False): # Generate contour plots for antagonistic parameters
	'''
	The primary method for producing contour plots. `retro` should in general be `DO_RETRO_CONTOUR` (it is only a parameter due to scoping
	issues).
	'''
	D.sel_advs = control_sa
	D.transm_probs = base_ta
	if ANTAG_TYPE not in antag_types: print('bad antagonistic type'); exit()
	antag_dims = antag_types[ANTAG_TYPE]
	fdir = f'a{ANTAG_TYPE}_contour'
	fulldir = f'full_outputs/{fdir}'
	param_to_change = D.sel_advs if ANTAG_TYPE == 'sel' else D.transm_probs
	antag_name = ANTAG_NMS[ANTAG_TYPE]
	h1_rng, vec_rng = [getRange(param_to_change[pop_nm], antag_dims[pop_nm], num_els=CONTOUR_DENSITY) for pop_nm in [H1, VEC]]
	lengths = {'Vector': 6, 'Host': 5}
	contour_data = {pop_nm: {(vec_idx, h1_idx): -1. for vec_idx, h1_idx in zip(vec_rng, h1_rng)} for pop_nm in lengths}
	round_idx = 0
	h1_bounds = [h1_rng[0], h1_rng[-1]]
	vec_bounds = [vec_rng[0], vec_rng[-1]]
	h1_pre_bounds = h1_bounds
	vec_pre_bounds = vec_bounds
	if DO_EXTEND:
		if not os.path.exists(fulldir): os.makedirs(fulldir)
		existing_dirs = [hv_dir for hv_dir in os.listdir(fulldir) if not ('png' in hv_dir or 'csv' in hv_dir) and '_' in hv_dir]
		h1s = []
		vecs = []
		for e_ds in existing_dirs:
			h1v, vecv = [float(hvv[1:]) for hvv in e_ds.split('_')]
			h1s += [h1v]
			vecs += [vecv]
		if existing_dirs:
			h1_pre_bounds = [min(h1s), max(h1s)]
			vec_pre_bounds = [min(vecs), max(vecs)]
		else:
			h1_pre_bounds = [0, 0]
			vec_pre_bounds = [0, 0]
	if not retro: # generate data if it's not already there
		for i in range(CONTOUR_DENSITY):
			for j in range(CONTOUR_DENSITY):
				round_idx += 1
				param_to_change[H1] = h1_rng[i]
				param_to_change[VEC] = vec_rng[j]
				fdir_ij = f'{fdir}/h{trunc(h1_rng[i])}_v{trunc(vec_rng[j])}'
				fulldir_ij = f'full_outputs/{fdir_ij}'
				num_inc_runs = NUM_RUNS - numSubDirs(fulldir_ij) if os.path.exists(fulldir_ij) else 0
				if num_inc_runs < 0: num_inc_runs = 0
				colorist.red(f'Beginning round {round_idx}/{CONTOUR_DENSITY**2} for antagonistic {antag_name}, {param_to_change}')
				if DO_EXTEND and ((isInRange(h1_rng[i], h1_pre_bounds) and isInRange(vec_rng[j], vec_pre_bounds)) and not num_inc_runs):
					colorist.red('\tSkipping round on account of extension')
					stat_lsts = {}
					if fdir_ij in os.listdir(fulldir):
						f = open(f'{fulldir_ij}/net_output_overall.opt', 'r')
						for opt_line in f.readlines():
							if 'I (' in opt_line:
								key, val = opt_line.split(':\t')
								stat_lsts[key] = [float(msv) for msv in val.split('\t+- ')]
				else: stat_lsts = doMultipleRuns(n=num_inc_runs, fdir=fdir_ij, force_multi=True, do_qc=True)
				for pop_nm, pop_len in lengths.items():
					if not stat_lsts: continue
					contour_data[pop_nm][(vec_rng[j], h1_rng[i])] = getNetAllFreq(stat_lsts, pop_len)
	for pop_nm, pop_data in contour_data.items():
		to_pop = [coord for coord, val in pop_data.items() if val < 0]
		for coord in to_pop: pop_data.pop(coord)
	colorist.cyan('Data generation complete. Reading data from full_outputs...')
	for opt_dir in os.listdir(fulldir):
		if ('png' in opt_dir or 'csv' in opt_dir) or '_' not in opt_dir: continue
		h1_c, vec_c = [float(hvc[1:]) for hvc in opt_dir.split('_')]
		if (h1_c, vec_c) in contour_data['Vector']: continue
		for pop_nm, pop_len in lengths.items():
			stat_lsts = doRetroStats(f'{fdir}/{opt_dir}')
			contour_data[pop_nm][(vec_c, h1_c)] = getNetAllFreq(stat_lsts, pop_len)
	colorist.cyan('Done.')
	for pop_nm, pop_data in contour_data.items():
		to_pop = [coord for coord in pop_data if not isInRange(coord[0], vec_bounds) or not isInRange(coord[1], h1_bounds)]
		[pop_data.pop(coord) for coord in to_pop]
	xtick_lbls = roundTrailing(*vec_rng)
	len_lbls = sum([len(str(xtl)) for xtl in xtick_lbls])
	if len_lbls > 50: xtick_lbls = [v for i, v in enumerate(xtick_lbls) if not i%2]
	xs_vec = []
	ys_h1 = []
	zs_all: dict[str, dict[tuple[float, float], float]] = {pop_nm: {} for pop_nm in FULL_NMS}
	
	def addLabels(colorbar_label: str='Mutated allele frequency'):
		param_type = ''
		if antag_name == 'selection': param_type = 'absolute fitness'
		elif antag_name == 'transmission': param_type = 'transmission probability'
		plt.colorbar(label=colorbar_label)
		plt.xlabel(f'Vector: {param_type}')
		plt.xticks(xtick_lbls, labels=xtick_lbls)
		plt.ylabel(f'Host: {param_type}')
	
	def plotScatter(pop_data: dict[tuple[float, float], float], save_fig: str, l_type: Literal['dec', 'pct']='dec', **kwargs):
		scatter_dir = f'{fulldir}/scatters'
		if not os.path.exists(scatter_dir): os.makedirs(scatter_dir)
		label_method = trunc
		def round_pct(f: float): return f'{roundNum(100*f)}%'
		if l_type == 'pct': label_method = round_pct
		plt.scatter(*transpose(pop_data.keys()), c=list(pop_data.values()), cmap=getColorMap(scale_type=COLOR_SCALE_TYPE))
		[plt.annotate(label_method(pop_data[coords]), coords, rotation=45) for coords in pop_data]
		addLabels(**kwargs)
		plt.savefig(f'{scatter_dir}/{save_fig}.png')
		plt.close()

	def getContourLine(pop_data: dict[tuple[float, float], float], mid_val: float=0.5, rng: float=0.05):
		ok_rng = getRange(mid_val, rng, num_els=2)
		coords_lst = [coords for coords, val in pop_data.items() if isInRange(val, ok_rng)]
		return coords_lst

	contour_kwargs = {'levels': 20, 'cmap': getColorMap(scale_type=COLOR_SCALE_TYPE)}
	lr_done = False
	for pop_nm in contour_data:
		fn_save = f'{fulldir}/{FULL_NMS[pop_nm]}'
		xs_vec, ys_h1, zs_freq = getGenContourData(contour_data[pop_nm], fn_save)
		if DO_EXTEND:
			nan_mask = np.isnan(zs_freq)
			xs_vec_mg, ys_h1_mg = np.meshgrid(xs_vec, ys_h1)
			xs_vec, ys_h1 = np.mgrid[min(xs_vec):max(xs_vec):100j, min(ys_h1):max(ys_h1):100j]
			zs_freq = np.array(zs_freq)
			zs_freq = griddata((xs_vec_mg[~nan_mask].flatten(), ys_h1_mg[~nan_mask].flatten()), zs_freq[~nan_mask].flatten(),
					  (xs_vec, ys_h1), method='cubic')
			for i in range(len(zs_freq)):
				for j in range(len(zs_freq[i])):
					if zs_freq[i,j] > 1: zs_freq[i,j] = 1
					elif zs_freq[i,j] < 0: zs_freq[i,j] = 0
		plt.contourf(xs_vec, ys_h1, zs_freq, **contour_kwargs)
		plt.title(f'{pop_nm}: {antag_name} contour plot ({SIM_LENGTH} days)')
		addLabels()
		plt.savefig(f'{fn_save}.png')
		if SHOW_CONTOUR: plt.show()
		plt.close()
		pop_data = contour_data[pop_nm]
		plotScatter(pop_data, save_fig=FULL_NMS[pop_nm])
		zs_all[pop_nm] = pop_data
		if ANTAG_TYPE == 'trans' and not lr_done:
			lr_base, lr_rng = 0.5, 0.05
			lr_coords_raw = getContourLine(pop_data, mid_val=lr_base, rng=lr_rng)
			lr_coords = transpose([(xc/vec_base_transm_p, yc/hst_base_transm_p) for (xc, yc) in lr_coords_raw])
			slope, intercept, r, p, stderr = linregress(*lr_coords)
			plt.scatter(*lr_coords)
			lr_xs = lr_coords[0]
			lr_ys = [slope*lrx + intercept for lrx in lr_xs]
			plt.plot(lr_xs, lr_ys, label=rf'slope -{roundAndSN(-slope)}, R$^2$ {roundAndSN(r**2)}, err {roundAndSN(stderr)}')
			plt.title(f'Linear regression: frequency {lr_base} (range {lr_rng})')
			plt.ylabel('Host: relative transmission advantage')
			plt.xlabel('Vector: relative transmission advantage')
			plt.legend()
			plt.savefig(f'{fn_save}_linreg.png')
			plt.show()
			plt.close()
			lr_done = True
	diff_data = {coord: zs_all['Host'][coord] - zs_all['Vector'][coord] for coord in zs_all['Host']}
	xv_diff, yh_diff, zs_diff = getGenContourData(diff_data, f'{fulldir}/diff')
	plt.contourf(xv_diff, yh_diff, zs_diff, **contour_kwargs)
	plt.title('Host-vector mutated allele frequency differences')
	addLabels(colorbar_label='Frequency difference')
	plt.savefig(f'{fulldir}/diff.png')
	plt.close()
	plotScatter(diff_data, save_fig='diff', l_type='pct')

def fixAntagDirs():
	'''
	Ensures the directories created by antagonistic parameters & associated contour plots are of the correct naming scheme. Originally
	used to switch from an older naming scheme to the current one; should not see much use overall.
	'''
	full_dir = f'full_outputs/a{ANTAG_TYPE}_contour'
	new_dir_names = {}
	test_str = f'{ANTAG_TYPE[0].upper()}{ANTAG_TYPE[1:]}'
	for hv_dir in os.listdir(full_dir):
		if '_' in hv_dir:
			opt_file = open(f'{full_dir}/{hv_dir}/net_output_overall.opt', 'r')
			for opt_line in opt_file.readlines():
				if test_str in opt_line:
					hv_nums = opt_line.split('{')[1][:-1]
					strs = hv_nums.split(', ')
					ndn = ''
					for str_spt in strs:
						pop_nm, pop_val = str_spt.split(': ')
						ndn += f'{pop_nm[1]}{trunc(pop_val)}_'
					new_dir_names[hv_dir] = ndn[:-2]
			opt_file.close()
	[os.rename(f'{full_dir}/{k}', f'{full_dir}/{v}') for k, v in new_dir_names.items()]

def getAllFreq(csv_data: dict[str, list[float]]) -> list[float]:
	div_fac = len(csv_data) - 3
	return list(np.sum(transpose([hdr.count(D.char)*np.array(vals)/div_fac for hdr, vals in csv_data.items() if 'I (' in hdr]), axis=1))

def compilePlots(output_dir: str, do_show: bool=True):
	data_dirs = getSubDirs(output_dir)
	times: list[float] = []
	freqs: dict[str, list[list[float]]] = {k: [] for k in FULL_NMS_R}
	for d_dir in data_dirs:
		full_dir = f'{output_dir}/{d_dir}'
		for pop_type in freqs:
			pop_data = readCSVData(f'{full_dir}/{pop_type}.csv')
			if not times: times = pop_data['times']
			freqs[pop_type] += [getAllFreq(pop_data)]
	for pop_type, pop_freqs in freqs.items():
		[plt.plot(times, run_freqs, label=i+1) for i, run_freqs in enumerate(pop_freqs)]
		plt.plot(times, 0*np.array(times), alpha=0.)
		plt.plot(times, [1 for t in times], alpha=0)
		plt.title(f'{FULL_NMS_R[pop_type]}: mutated allele frequency ({len(pop_freqs)} simulations)')
		plt.xlabel('Simulation time (days)')
		plt.ylabel('Mutated allele frequency')
		plt.legend()
		plt.savefig(f'{output_dir}/{pop_type}.png')
		if do_show: plt.show()
		plt.close()

if __name__ == '__main__':
	# Change the below string to pick which case to run

	match TO_RUN:
		### For computational cost tracking
		case 'cost':
			doTimeBreakdown()

		### Primary
		case 'primary':
			if DO_ANTAG: doContourPlots(retro=DO_RETRO_CONTOUR); exit()
			fdir = 'temp_dir'
			if FILE_DIR: fdir = FILE_DIR
			num_comp_runs = numSubDirs(f'full_outputs/{fdir}')
			doMultipleRuns(n=NUM_RUNS-num_comp_runs, fdir=fdir, force_multi=(NUM_RUNS != 1))

		### To get overall data from the given pre-generated directory
		case 'overall':
			doRetroStats(FILE_DIR)
			compilePlots(f'full_outputs/{FILE_DIR}')