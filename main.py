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
    'do_mixed_infs': False,
}
HST1 = {
    'bd': 0.,
    'ir': 0.,
    'rr': .05,
    'wi': .20,
    'pn': 'h1',
    'pn_full': 'Host 1',
    'do_mixed_infs': False,
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
    'pc': 12,
    'mut_chance': 4e-2,
    'para_lsp': 2.,
    'is_hap': False,
    'do_sr': False,
    'do_mutation': False,
    'do_indvs': True,
    'pc_to_transmit': 6,
}
INDV_HST = {
    'pc': 120,
    'mut_chance': 8e-5,
    'para_lsp': 2.,
    'is_hap': True,
    'do_sr': False,
    'do_mutation': True,
    'do_indvs': False,
    'pc_to_transmit': 60,
}
INDVS = [INDV_VEC, INDV_HST]

# for p_fac of 5e4, nt 2e4: epidemic params are 4e3, 1e3, 7e1 for ir, rr, wi respectively (stab/epi)

A = allele(char='A', fav_pop='h1', unf_pop='h2', param='itr', fac=0.3)
B = allele(char='B', fav_pop='h1', unf_pop='h2', param='itr', fac=0.3)
C = allele(char='C', fav_pop='h1', unf_pop='h2', param='rr', fac=-0.6)
D = allele(char='D', fav_pop='h1', unf_pop='h2', param='itr', fac=1.25)
ALLELES = [D]

PARAMS_1 = HST1
PARAMS_2 = VEC
PARAMS_3 = HST2

def run(p0: np.ndarray=np.array([[20, 1, 0], [21, 0, 0]], dtype='float64'), p_fac: float=100., nt: float=2.,
        plot_res: bool=True, t_scale: float=200., do_allele_mod: bool=True, weight_infs: bool=True,):
    '''
    Run the simulation.

    ### Parameters
    - `p0`: The initial population ratios (S, I, R) as a NumPy array of 3-element NumPy arrays.
    - `p_fac`: The scale factor on population.
    - `nt`: Time steps per day. Alternatively, 24/`nt` is hours per time step.
    - `plot_res`: Whether or not to display the results' graph, as a bool.
    - `t_scale`: The number of days to run the simulation for.
    - `do_allele_mod`: Whether or not to use the allele-based mutation model, as a bool.
    - `is_haploid`: Whether the model is haploid (AbC) or diploid (AabbCC).
    '''
    [os.remove(file) for file in os.listdir() if file.endswith('.dat')]
    f = open('inf_events_raw.dat', 'x')
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
        i_params['store_chance'] = 1e1/(p_fac*t_scale*nt*i_params['pc'])
    nt = float(int(nt*t_scale))
    # [p0_1, p0_2, p0_3] = [population(p0[i], **INDV) for i in range(3)]
    hosts_1 = population(p0[0], **INDV_HST)
    vectors = population(p0[1], **INDV_VEC)
    m1 = SIR(hosts_1, **PARAMS_1)
    m2 = SIR(vectors, **PARAMS_2)
    # m3 = SIR(p0_3, **PARAMS_3)
    itr_h1 = 0.15 # h1 <-> vec
    itr_h2 = 0.10 # h2 <-> vec
    # m1.itr = {p0_2: itr_h1, p0_3: 0.} # the number represents the rate at which m1 infects that population
    # m2.itr = {p0_1: itr_h1, p0_3: itr_h2}
    # m3.itr = {p0_1: 0., p0_2: itr_h2}
    m1.itr = {vectors: itr_h1}
    m2.itr = {hosts_1: itr_h1}
    t0 = time.time()
    # mdls = [m1, m2, m3]
    mdls = [m1, m2]
    ts, ps, times, pops, ps_unwgt = simShell(t_max, mdls, nt=nt, alleles=alleles, weight_infs=weight_infs)
    ex_tm = time.time() - t0
    times_norm = list(100*normalise(np.array(times)))
    print(f'Execution time: {ex_tm}')
    print('Breakdown:')
    [print(f'{i}:\t{times_norm[i]}') for i in range(len(times))]
    print(f'Extra time: {ex_tm - sum(times)}')
    [p.printDat() for p in pops]
    # dimensions of ps: layer 1 is times, layer 2 is models at that time, layer 3 is pop #s for that model
    settings = {}
    settings['hap/dip'] = 'diploid'
    if INDV_HST['is_hap']: settings['hap/dip'] = 'haploid'
    if INDV_HST['is_hap'] and not INDV_VEC['is_hap']: settings['hap/dip'] = 'hap./dip. hybrid'
    settings['sr'] = 'sexual reproduction'
    if not INDV_VEC['do_sr']: settings['sr'] = f'no {settings["sr"]}'
    settings['mut'] = 'no'
    if INDV_HST['do_mutation']: settings['mut'] = 'host'
    settings['mut'] += ' mutation'
    f = open('inf_events.dat', 'x')
    f_raw = open('inf_events_raw.dat', 'r')
    inf_events = {}
    for line in f_raw.readlines():
        if line in inf_events: inf_events[line] += 1
        else: inf_events[line] = 1
    f.write(str(inf_events))
    f.close()
    f_raw.close()
    
    plots = {} # title (str) : plt
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
    # getDims(ps)
    for i in range(len(mdls)):
        ns = [''.join(n.split('.')) for n in pops[i].getAllPopNms()]
        gens = []
        for n in ns:
            if '(' in n and ')' in n: gens += [n[n.index('(')+1:n.index(')')]]
            else: gens += [n]
        ps_i = np.array([k[i] for k in ps])
        [plt.plot(ts, ps_i[:,j], label=ns[j], color=str2Color(gens[j]), alpha=pop2Alpha(ns[j])) for j in range(len(ns))
         if (ns[j][0] != 'R') and (ns[j][0] != 'S')]
        ps_unwgt_i = np.array([k[i] for k in ps_unwgt])
        if weight_infs and mdls[i].pop.do_indvs: [plt.plot(ts, ps_unwgt_i[:,j], label=f'{ns[j]} (unweighted)', color=str2Color(gens[j]),
            alpha=pop2Alpha(ns[j])/2) for j in range(len(ns)) if (ns[j][0] != 'R') and (ns[j][0] != 'S')]
        # plt.plot(ts, sum(ps_i.transpose()), label='N')
        plt.title(f'{mdls[i].pn_full} population ({", ".join(list(settings.values()))})')
        plt.legend()
        plt.xlabel('Simulation time')
        plt.ylabel('Population')
        plt.savefig(f'{fn(mdls[i].pn)}.png')
        fig = plt.gcf()
        plots[fig.axes[0].get_title()] = fig
        if plot_res: plt.show()
        plt.close()
    return plots

if __name__ == '__main__':
    run()
    # plots: dict[str, plt.Figure] = {}
    # INDV['do_sr'] = True
    # plots.update(run())
    # INDV['do_sr'] = False
    # plots.update(run())
    # plots_new: dict[str, list[plt.Figure]] = {}
    # for p in plots:
    #     pop = p.split(' ')[0]
    #     if pop in plots_new: plots_new[pop] += [plots[p]]
    #     else: plots_new[pop] = [plots[p]]
    # for p in plots_new:
    #     for f in plots_new[p]:
    #         [plt.plot(ax.lines[0].get_xdata(), ax.lines[0].get_ydata(), label=f'{ax.get_label()} ({ax.get_title()})') for ax in f.axes]
    #         plt.title(f'Combined graph (pop {p})')
    #         plt.legend()
    #         plt.xlabel(f.axes[0].get_xlabel())
    #         plt.ylabel(f.axes[0].get_ylabel())
    #         plt.savefig(f'comb_{p}.png')
    #         plt.show()
    #         plt.close()