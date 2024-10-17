from sim_lib import *
from allele import *
from color import *
from time import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import cProfile
import pstats

# Rates are in terms of '# events expected per day'
VEC = {             # Parameters for transmission vector behavior (mosquito)
    'bd': 0.071,    # Birth/death rate
    'ir': 0.,       # Infection rate between individuals. Largely obsolete for a vector-based model
    'rr': 0,        # Recovery rate
    'wi': 0.,       # Waning immunity rate
    'pn': 'vec',            # Short population name, used in console output & within the code
    'pn_full': 'Vector',    # Full population name, used in graph titles
    'is_vector': True,      # Whether or not the population is a disease vector
}

tau = 90                    # Average duration of immunity (days)
hst_base_transm_p = 0.11    # Probability of contact resulting in a host-vector transmission
ct_rate = 0.28              # Contact rate (average bites per day)
h_ir = ct_rate*hst_base_transm_p                            # Infection rate (only used in the context of wi_rate)
wi_rate = h_ir*np.exp(-h_ir*tau)/(1 - np.exp(-h_ir*tau))    # Waning immunity rate

HST1 = {
    'bd': 0.,
    'ir': 0.,
    'rr': 2.19e-3,
    'wi': wi_rate,
    'pn': 'h1',
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

all_adv = 1.5      # Selection advantage parameter
wld_adv = 1/all_adv # Pro-wild allele (lowercase) selection advantage. Represented as a disadvantage for the mutated allele
mut_adv = all_adv   # Pro-mutated allele (uppercase) selection advantage
D.sel_advs = {'h1': 1.0, 'vec': mut_adv} # Selection biases for mutated allele, relative to the wild allele (always 1.0)

NUM_RUNS = 5                        # Number of simulations to run
FILE_DIR = 'seladv_vec_explicit'    # Leave blank for a procedural directory name

# For transmission probabilities: the pop ID is the source -- i.e., 'vec': 0.021 means 2.1% transmission chance from vector to host
D.base_transm_probs = {'h1': hst_base_transm_p, 'vec': 0.021} # for wild-type allele

# Switch between the two lines below to give transmission advantages/disadvantages
D.transm_probs = D.base_transm_probs.copy()
# D.transm_probs = {'h1': hst_base_transm_p, 'vec': 0.07}

ALLELES = [D] # Do NOT have more than one allele here -- the simulation has been optimised for the single-locus case.
              # Adding more WILL break it!

PARAMS_1 = HST1
PARAMS_2 = VEC

def run(p0: np.ndarray=np.array([[20, 1, 0], [21, 0, 0]], dtype='float64'), p_fac: float=1200., nt: float=1., num_hist: int=0,
        plot_res: bool=True, t_scale: float=20000., init_mut_prop: float=0.5, fdir: str=''):
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
    mkFile('inf_events_raw.dat', 'last_event.dat',)
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
    
    frac_to_take = 0.2                  # Fraction of the simulation to use when recording means (in net_output.opt)
    output_fn = f'{fdir}net_output.opt' # '.opt' is a plain text file. Only marked that way to keep the file clearing from catching it
    if os.path.exists(output_fn): os.remove(output_fn)
    f = open(output_fn, 'x')
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
                data_to_keep = plt_lst[int(-frac_to_take*len(plt_lst)):]
                f.write(f'{ns[j]}:\t{roundAndSN(np.mean(data_to_keep))}\t+- {roundAndSN(np.std(data_to_keep))}\n')
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
    hists_with_pn = {'vec': hists_v, 'h1': hists_h}
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
            if hpn == 'h1': plt.ylim(top=hist_max)
            plt.savefig(f'hists/{hpn}_{i}.png')
            plt.close()

def doTimeBreakdown(): # Runs the simulation with a profiler & processes the output accordingly
    time_output_fn = 'time_breakdown.dat'
    cProfile.run('run()', time_output_fn)
    time_output_txt = 'time_breakdown.txt'
    time_output_file = open(time_output_txt, 'x')
    p_stats = pstats.Stats(time_output_fn, stream=time_output_file)
    p_stats.strip_dirs()
    p_stats.sort_stats('cumtime')
    p_stats.print_stats()

def doMultipleRuns(n: int=3, fdir: str=''): # Run the simulation multiple times & save the results to a particular directory
    if n == 1: run(); return
    full_dir = f'full_outputs/{fdir}/'
    for i in range(1,n+1):
        print(10*'-' + f' RUN {i} ' + 10*'-')
        run(fdir=f'{full_dir}{i}/')

if __name__ == '__main__':
    # run()
    # doTimeBreakdown()

    adv_amt = str(mut_adv).split('.')[-1]
    adv_tgt = 'host'
    adv_type = 'sel'
    fdir = f'{adv_type}adv_{adv_tgt}_{adv_amt}'
    if FILE_DIR: fdir = FILE_DIR
    doMultipleRuns(n=NUM_RUNS, fdir=fdir)