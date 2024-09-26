from sim_lib import *
from allele import *
from color import *
import matplotlib.pyplot as plt
import time
import os

### Deprecated parameters
STAB = {        # Parameters for stable behavior (equilibrium)
    'bd': 0,        # Birth/death rate
    'ir': 60,       # Infection rate
    'rr': 28,       # Recovery rate
    'wi': 10.5,     # Waning immunity rate
    'nm': 'Stable'
}
EPI = {         # Parameters for epidemic behavior (short-lived spikes of infection)
    'bd': 0,        # Birth/death rate
    'ir': 4e3,      # Infection rate
    'rr': 1e3,      # Recovery rate
    'wi': 7e1,       # Waning immunity rate
    'nm': 'Epidemic'
}

### Up-to-date parameters
# Note: The interspecific transmission rate is in the main method itself!
# Rates are in terms of '# events expected per day'
VEC = {         # Parameters for transmission vector behavior (mosquito)
    'bd': 0.3,  # Birth/death rate
    'ir': 0.,   # Infection rate between individuals. Largely obsolete for a vector-based model
    'rr': 0,    # Recovery rate
    'wi': 0.,   # Waning immunity rate
    'pn': 'vec',            # Short population name, used in console output & within the code
    'pn_full': 'Vector',    # Full population name, used in graph titles
    'is_vector': True,      # Whether or not the population is a disease vector
    'do_mixed_infs': True,
}
HST1 = {
    'bd': 0.,
    'ir': 0.,
    'rr': .05,
    'wi': .20,
    'pn': 'h1',
    'pn_full': 'Host 1',
    'do_mixed_infs': True,
}
HST2 = {
    'bd': 0.,
    'ir': 0.,
    'rr': .05,
    'wi': .20,
    'pn': 'h2',
    'pn_full': 'Host 2',
}

INDV_VEC = {
    'pc': 120,
    'mut_chance': 4e-3,
    'para_lsp': 2.,
    'is_hap': False,
    'do_sr': False,
    'do_mutation': False,
    'do_indvs': True,
    'do_sel_bias': True,
}
INDV_HST = {
    'pc': int(1.2e6),
    'mut_chance': 8e-6,
    'para_lsp': 2.,
    'is_hap': True,
    'do_sr': False,
    'do_mutation': True,
    'do_indvs': True,
    'do_sel_bias': True,
}
INDVS = [INDV_VEC, INDV_HST]
for INDV in INDVS: INDV['pc_to_transmit'] = int(INDV['pc']/2)

# for p_fac of 5e4, nt 2e4: epidemic params are 4e3, 1e3, 7e1 for ir, rr, wi respectively (stab/epi)

A = allele(char='A', fav_pop='h1', unf_pop='h2', param='itr', fac=0.3)
B = allele(char='B', fav_pop='h1', unf_pop='h2', param='itr', fac=0.3)
C = allele(char='C', fav_pop='h1', unf_pop='h2', param='rr', fac=-0.6)
D = allele(char='D', fav_pop='h1', unf_pop='h2', param='itr', fac=0.0)

mut_adv = 1.2
wld_adv = 1/mut_adv
D.sel_advs = {'h1': 1.05, 'vec': 1.0}
D.trans_advs = {'h1': 1.0, 'vec': 1.0}

ALLELES = [D]

PARAMS_1 = HST1
PARAMS_2 = VEC
PARAMS_3 = HST2

def run(p0: np.ndarray=np.array([[20, 1, 0], [21, 0, 0]], dtype='float64'), p_fac: float=1200., nt: float=2., num_hist: int=0,
        plot_res: bool=False, t_scale: float=100., do_allele_mod: bool=True, weight_infs: bool=True, do_mix_start: bool=False,):
    '''
    Run the simulation.

    ### Parameters
    - `p0`: The initial population ratios (S, I, R) as a NumPy array of 3-element NumPy arrays.
    - `p_fac`: The scale factor on population.
    - `nt`: Time steps per day. Alternatively, 24/`nt` is hours per time step.
    - `plot_res`: Whether or not to display the results' graph, as a bool.
    - `t_scale`: The number of days to run the simulation for.
    - `do_allele_mod`: Whether or not to use the allele-based mutation model, as a bool.
    - `weight_infs`: Whether or not to display 'weighted' data for infections (weighted according to the strains' relative
        genotype frequencies).
    - `do_mix_start`: Whether or not to have a mixed distribution of infected individuals (wild & mutated) or uniform (just wild).
    '''
    [os.remove(file) for file in os.listdir() if file.endswith('.dat')]
    if 'hists' not in os.listdir(): os.mkdir('hists')
    f = open('inf_events_raw.dat', 'x')
    f.close()
    f = open('last_event.dat', 'x')
    f.close()
    alleles = []
    if do_allele_mod: alleles = ALLELES
    p0 *= p_fac
    t_max = t_scale
    for i in range(len(INDVS)):
        i_params = INDVS[i]
        para_gens = int((24/nt)/i_params['para_lsp'])
        i_params['para_gens'] = para_gens
        i_params['alleles'] = alleles
        # i_params['store_chance'] = 5e1/(p_fac*t_scale*nt*i_params['pc'])
    nt = float(int(nt*t_scale))
    # [p0_1, p0_2, p0_3] = [population(p0[i], **INDV) for i in range(3)]
    hosts_1 = population(p0[0], **INDV_HST)
    vectors = population(p0[1], **INDV_VEC)
    m1 = SIR(hosts_1, **PARAMS_1)
    m2 = SIR(vectors, **PARAMS_2)
    # m3 = SIR(p0_3, **PARAMS_3)
    itr_h1 = 0.25 # h1 <-> vec
    itr_h2 = 0.10 # h2 <-> vec
    # m1.itr = {p0_2: itr_h1, p0_3: 0.} # the number represents the rate at which m1 infects that population
    # m2.itr = {p0_1: itr_h1, p0_3: itr_h2}
    # m3.itr = {p0_1: 0., p0_2: itr_h2}
    m1.itr = {vectors: itr_h1}
    m2.itr = {hosts_1: itr_h1}
    t0 = time.time()
    # mdls = [m1, m2, m3]
    mdls = [m1, m2]
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
    
    def getDims(lst: list, tab: str=''):
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
    
    # print(hists_h)
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