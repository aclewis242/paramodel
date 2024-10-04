from sim_lib import *
from allele import *
from color import *
import matplotlib.pyplot as plt
import time
import os

# Rates are in terms of '# events expected per day'
VEC = {         # Parameters for transmission vector behavior (mosquito)
    'bd': 0.3,  # Birth/death rate
    'ir': 0.,   # Infection rate between individuals. Largely obsolete for a vector-based model
    'rr': 0,    # Recovery rate
    'wi': 0.,   # Waning immunity rate
    'pn': 'vec',            # Short population name, used in console output & within the code
    'pn_full': 'Vector',    # Full population name, used in graph titles
    'is_vector': True,      # Whether or not the population is a disease vector
}
HST1 = {
    'bd': 0.,
    'ir': 0.,
    'rr': .05,
    'wi': .20,
    'pn': 'h1',
    'pn_full': 'Host 1',
}
BTP = 0.5
ITR = 0.25 *(1/BTP)

INDV_VEC = {
    'pc': 120,
    'mut_chance': 4e-3,
    'para_lsp': 2.,
    'is_hap': False,
    'do_mutation': False,
}
INDV_HST = {
    'pc': int(1.2e6),
    'mut_chance': 8e-6,
    'para_lsp': 2.,
    'is_hap': True,
    'do_mutation': True,
}
INDVS = [INDV_VEC, INDV_HST]
for INDV in INDVS: INDV['pc_to_transmit'] = int(INDV['pc']/2)

D = allele(char='D')

mut_adv = 1.05
wld_adv = 1/mut_adv
D.sel_advs = {'h1': mut_adv, 'vec': 1.0}
transm_p = 0.5
D.transm_probs = {'h1': transm_p, 'vec': transm_p} # pop ID is the source -- e.g. 'h1' means 'prob of transmission from h1'
D.base_transm_probs = {'h1': BTP, 'vec': BTP} # for wild-type allele

ALLELES = [D] # Do NOT have more than one allele here -- the simulation has been optimised for the single-locus case.
              # Adding more WILL break it!

PARAMS_1 = HST1
PARAMS_2 = VEC

def run(p0: np.ndarray=np.array([[20, 1, 0], [21, 0, 0]], dtype='float64'), p_fac: float=1200., nt: float=2., num_hist: int=0,
        plot_res: bool=False, t_scale: float=100., weight_infs: bool=True, do_mix_start: bool=False,):
    '''
    Run the simulation.

    ### Parameters
    - `p0`: The initial population ratios (S, I, R) as a NumPy array of 3-element NumPy arrays.
    - `p_fac`: The scale factor on population.
    - `nt`: Time steps per day. Alternatively, 24/`nt` is hours per time step.
    - `num_hist`: The number of histograms to generate (varying over time).
    - `plot_res`: Whether or not to display the results' graph, as a bool.
    - `t_scale`: The number of days to run the simulation for.
    - `weight_infs`: Whether or not to display 'weighted' data for infections (weighted according to the strains' relative
        genotype frequencies).
    - `do_mix_start`: Whether or not to have a mixed distribution of infected individuals (wild & mutated) or uniform (just wild).
    '''
    [os.remove(file) for file in os.listdir() if file.endswith('.dat')]
    mkDir('hists', 'old images')
    mkFile('inf_events_raw.dat', 'last_event.dat',)
    alleles = ALLELES
    p0 *= p_fac
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
    m1.itr = {vectors: ITR}
    m2.itr = {hosts_1: ITR}
    mdls = [m1, m2]
    t0 = time.time()
    ts, ps, times, pops, ps_unwgt, vpis, hpis, hists_v, hists_h, hist_tms = simShell(
        t_max, mdls, nt=nt, alleles=alleles, weight_infs=weight_infs, do_mix_start=do_mix_start, num_hist=num_hist)
    ex_tm = time.time() - t0
    times_norm = normPercentList(times)
    print(f'\nExecution time: {roundNum(ex_tm, prec=3)}') # consider colored console output for readability
    print('Breakdown:')
    printFloatList(times_norm)
    print(f'Extra time: {ex_tm - sum(times)}')
    # for p in pops:
    #     print(f'{p.pn} time breakdown:')
    #     printFloatList(normPercentList(p.times))
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
    
    def getDims(lst: list, tab: str=''): # Displays the dimensions of the given list. Useful when handling complex/unorthodox structures.
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
    for i in range(len(mdls)):
        ns = [''.join(n.split('.')) for n in pops[i].getAllPopNms()]
        gens = []
        for n in ns:
            if '(' in n and ')' in n: gens += [n[n.index('(')+1:n.index(')')]]
            else: gens += [n]
        # ps_unwgt_i = np.array([k[i] for k in ps_unwgt])
        # if weight_infs and mdls[i].pop.do_indvs: [plt.plot(ts, ps_unwgt_i[:,j], label=f'{ns[j]} (unweighted)', color=str2Color(gens[j]),
        #     alpha=pop2Alpha(ns[j])/4) for j in range(len(ns)) if (ns[j][0] != 'R') and (ns[j][0] != 'S')]
        ps_i = np.array([k[i] for k in ps])
        [plt.plot(ts, ps_i[:,j], label=ns[j], color=str2Color(gens[j]), alpha=pop2Alpha(ns[j])) for j in range(len(ns))
         if (ns[j][0] != 'R') and (ns[j][0] != 'S')]
        net_i = vpis
        if mdls[i].pn == 'h1': net_i = hpis
        plt.plot(ts, net_i, label='I (total)')
        plt.plot(ts, 0*ts, alpha=0.)
        plt.title(f'{mdls[i].pn_full} population')
        plt.legend()
        plt.xlabel('Simulation time')
        plt.ylabel('Population')
        plt.savefig(f'{fn(mdls[i].pn)}.png')
        if plot_res: plt.show()
        plt.close()
    
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

if __name__ == '__main__':
    run()