from sim_lib import *
from allele import *
from color import *
from time import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import cProfile
import pstats

H1 = 'h1'
VEC = 'vec'

# Rates are in terms of '# events expected per day'
VCT = {             # Parameters for transmission vector behavior (mosquito)
    'bd': 0.071,    # Birth/death rate
    'ir': 0.,       # Infection rate between individuals. Largely obsolete for a vector-based model
    'rr': 0,        # Recovery rate
    'wi': 0.,       # Waning immunity rate
    'pn': VEC,              # Short population name, used in console output & within the code
    'pn_full': 'Vector',    # Full population name, used in graph titles
    'is_vector': True,      # Whether or not the population is a disease vector
}

hst_transm_p_wi = 0.11      # The base transmission probability for hosts, only used here to come up with wi_rate
ct_rate = 0.28              # Contact rate (average bites per day)
tau = 90                    # Average duration of immunity (days)
h_ir = ct_rate*hst_transm_p_wi                              # Infection rate (only used in the context of wi_rate)
wi_rate = h_ir*np.exp(-h_ir*tau)/(1 - np.exp(-h_ir*tau))    # Waning immunity rate

HST1 = {
    'bd': 0.,
    'ir': 0.,
    'rr': 2.19e-3,
    'wi': wi_rate,
    'pn': H1,
    'pn_full': 'Host',
}

para_lifespan = 8.              # Parasite lifespan (hours)
INDV_VEC = {
    'pc': 20,                   # Parasite count
    'para_lsp': para_lifespan, 
    'is_hap': False,            # Whether or not this individual's parasites are haploid
    'do_mutation': False,       # Whether or not this individual's parasites can mutate
    'pc_to_transmit': 15,       # The number of parasites transmitted per infection
}
INDV_HST = {
    'pc': int(1e8),
    'mut_chance': 2.94e-6,      # Chance of mutation per parasite per generation
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
control_sa = {H1: 1.0, VEC: 1.0}        # Used for the control case
hst_sa_inv = {H1: 1.05, VEC: 1.0}       # Used for the invasion example of the host selection advantage case
hst_sa_exp = {H1: 1.005, VEC: 1.0}      # Used for the host selection advantage case (quantitative)
vec_sa_exp = {H1: 1.0, VEC: 1.5}        # Used for the vector selection advantage case
hv_sel_ant = {H1: 1.005, VEC: 1/1.5}    # Used for the host-vector antagonistic selection case
vh_sel_ant = {H1: 1/1.005, VEC: 1.5}    # Used for the vector-host antagonistic selection case

D.sel_advs = control_sa

NUM_RUNS = 3            # Number of simulations to run
FILE_DIR = 'test_opt'   # Directory to save files under, if multiple simulations are being run
INIT_MUT_PROP = 0.5     # Initial proportion of mutated alleles
SIM_LENGTH = 500       # The length of the simulation (days)
SHOW_RES = True         # Whether or not to show the results
PROP_FOR_STAT = 0.2     # The proportion (from the end) of the results to use for statistics

# For transmission probabilities: the pop ID is the source -- i.e., vec: 0.021 means 2.1% transmission chance from vector to host

hst_base_transm_p = 0.11    # Probability of contact resulting in a host-vector transmission
hst_adv_transm_p = 0.12     # Advantageous transmission probability for hosts
hst_dis_transm_p = 0.10     # Disadvantageous transmission probability for hosts

vec_base_transm_p = 0.021   # Probability of contact resulting in a vector-host transmission
vec_adv_transm_p = 0.023    # Advantageous transmission probability for vectors
vec_dis_transm_p = 0.019    # Disadvantageous transmission probability for vectors

D.base_transm_probs = {H1: hst_base_transm_p, VEC: vec_base_transm_p} # for wild-type allele

base_ta = D.base_transm_probs.copy()                    # No advantages
hst_ta = {H1: hst_adv_transm_p, VEC: vec_base_transm_p} # Host transmission advantage
vec_ta = {H1: hst_base_transm_p, VEC: vec_adv_transm_p} # Vector transmission advantage
hv_at = {H1: hst_adv_transm_p, VEC: vec_dis_transm_p}   # Host-vector antagonistic transmission
vh_at = {H1: hst_dis_transm_p, VEC: vec_adv_transm_p}   # Vector-host antagonistic transmission

D.transm_probs = hst_ta

##################################
### ----- END USER INPUT ----- ###
##################################

ALLELES = [D] # Do NOT have more than one allele here -- the simulation has been optimised for the single-locus case.
              # Adding more WILL break it!

PARAMS_1 = HST1
PARAMS_2 = VCT

def run(p0: np.ndarray=np.array([[20, 1, 0], [21, 0, 0]], dtype='float64'), p_fac: float=1200., nt: float=1., num_hist: int=0,
        plot_res: bool=SHOW_RES, t_scale: float=SIM_LENGTH, init_mut_prop: float=INIT_MUT_PROP, fdir: str=''):
    '''
    Run the simulation.

    ### Parameters
    - `p0`: The initial population ratios (S, I, R) as a NumPy array of 3-element NumPy arrays.
    - `p_fac`: The scale factor on population.
    - `nt`: Time steps per day. Alternatively, 24/`nt` is hours per time step.
    - `num_hist`: The number of histograms to generate (varying over time).
    - `plot_res`: Whether or not to display the results' graph, as a bool.
    - `t_scale`: The number of days to run the simulation for.
    - `init_mut_prop`: The initial fraction of mutated (uppercase) alleles. Must be between 0 and 1.
    - `fdir`: The directory to save the results under.
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
    ts, ps, pops, hists_v, hists_h, hist_tms = simShell(t_max, mdls, nt=nt, alleles=alleles, init_mut_prop=init_mut_prop, num_hist=num_hist)
    ex_tm = time() - t0
    print(f'\nExecution time: {roundNum(ex_tm, prec=3)}') # consider colored console output for readability

    [p.printDat() for p in pops]
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
                mean_stdev = saveStats(plt_lst, frac_to_take=PROP_FOR_STAT)
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
        plt.savefig(f'{file_nm}_freq.png')
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
    # structure of all_data: row nm (I (..)) : [[mean 1, std 1], [mean 2, std 2], ...]
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
    f.close()

def doMultipleRuns(n: int=3, fdir: str=''): # Run the simulation multiple times & save the results to a particular directory
    if n == 1: run(); return
    full_dir = f'full_outputs/{fdir}/'
    all_data: dict[str, list[list[float]]] = {}
    for i in range(1,n+1):
        print(10*'-' + f' RUN {i} ' + 10*'-')
        ret_data = run(fdir=f'{full_dir}{i}/')
        for k in ret_data:
            if k in all_data: all_data[k] += [ret_data[k]]
            else: all_data[k] = [ret_data[k]]
    doOverallData(all_data, full_dir)

def doRetroStats(fdir: str):
    full_dir = f'full_outputs/{fdir}/'
    csv_names = {'Vector': VEC, 'Host': H1}
    bad_cols = ['Unnamed', 'times']
    all_data: dict[str, list[list[float]]] = {}
    for r_dir in os.listdir(full_dir):
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
                stats_to_add = list(saveStats(v, frac_to_take=PROP_FOR_STAT))
                if k in all_data: all_data[k] += [stats_to_add]
                else: all_data[k] = [stats_to_add]
    doOverallData(all_data, full_dir)

if __name__ == '__main__':
    # run()
    # doTimeBreakdown()

    fdir = 'temp_dir'
    if FILE_DIR: fdir = FILE_DIR
    doMultipleRuns(n=NUM_RUNS, fdir=fdir)

    # doRetroStats('seladv_vec_invasion')