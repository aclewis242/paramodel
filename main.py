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
    'mr': 1e-2,             # Mutation rate
    'pn': 'vec',            # Short population name, used in console output & within the code
    'pn_full': 'Vector',    # Full population name, used in graph titles
    'is_vector': True,      # Whether or not the population is a disease vector
}
HST1 = {
    'bd': 0.,
    'ir': 0.,
    'rr': .05,
    'wi': .20,
    'mr': 1e-3,
    'pn': 'h1',
    'pn_full': 'Host 1',
}
HST2 = {
    'bd': 0.,
    'ir': 0.,
    'rr': .05,
    'wi': .20,
    'mr': 1e-3,
    'pn': 'h2',
    'pn_full': 'Host 2',
}

# for p_fac of 5e4, nt 2e4: epidemic params are 4e3, 1e3, 7e1 for ir, rr, wi respectively (stab/epi)

A = allele(char='A', fav_pop='h2', unf_pop='h1', param='itr', fac=2.0)
B = allele(char='B', fav_pop='h1', unf_pop='h2', param='wi', fac=0.5)
C = allele(char='C', fav_pop='h1', unf_pop='h2', param='rr', fac=-0.6)
ALLELES = [A]

PARAMS_1 = HST1
PARAMS_2 = VEC
PARAMS_3 = HST2

def run(p0: np.ndarray=np.array([[20, 1, 0], [21, 0, 0], [20, 1, 0]], dtype='float64'), p_fac: float=50, nt: float=6.,
        plot_res: bool=True, t_scale: float=500., do_allele_mod: bool=True, is_haploid: bool=True, para_lsp: float=2., pc: int=12,
        is_hyb: bool=True):
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
    alleles = []
    if do_allele_mod: alleles = ALLELES
    p0 *= p_fac
    t_max = t_scale
    para_gens = int((24/nt)/para_lsp)
    nt = float(int(nt*t_scale))
    [p0_1, p0_2, p0_3] = [population(p0[i], pc=pc, i_h=is_haploid) for i in range(3)]
    m1 = SIR(p0_1, **PARAMS_1)
    m2 = SIR(p0_2, **PARAMS_2)
    m3 = SIR(p0_3, **PARAMS_3)
    itr_h1 = 0.30 # h1 <-> vec
    itr_h2 = 0.10 # h2 <-> vec
    m1.itr = {p0_2: itr_h1, p0_3: 0.} # the number represents the rate at which m1 infects that population
    m2.itr = {p0_1: itr_h1, p0_3: itr_h2}
    m3.itr = {p0_1: 0., p0_2: itr_h2}
    t0 = time.time()
    mdls = [m1, m2, m3]
    if is_hyb:
        for m in mdls:
            if m.is_vector: m.pop.is_dip = True
    ts, ps, times, pops = simShell(t_max, mdls, nt, alleles, para_gens)
    ex_tm = time.time() - t0
    times_norm = list(100*normalise(np.array(times)))
    print(f'Execution time: {ex_tm}')
    print('Breakdown:')
    [print(f'{i}:\t{times_norm[i]}') for i in range(len(times))]
    print(f'Extra time: {ex_tm - sum(times)}')
    [p.printDat() for p in pops]
    # dimensions of ps: layer 1 is times, layer 2 is models at that time, layer 3 is pop #s for that model
    for i in range(len(mdls)):
        ns = [''.join(n.split('.')) for n in pops[i].getAllPopNms()]
        gens = []
        for n in ns:
            if '(' in n and ')' in n: gens += [n[n.index('(')+1:n.index(')')]]
            else: gens += [n]
        ps_i = np.array([k[i] for k in ps])
        [plt.plot(ts, ps_i[:,j], label=ns[j], color=str2Color(gens[j]), alpha=pop2Alpha(ns[j])) for j in range(len(ns)) if ns[j][0] != 'R']
        plt.plot(ts, sum(ps_i.transpose()), label='N')
        plt.title(f'{mdls[i].pn_full} population')
        plt.legend()
        plt.xlabel('Simulation time')
        plt.ylabel('Population')
        plt.savefig(f'{fn(mdls[i].pn)}.png')
        if plot_res: plt.show()
        plt.close()
    return ps

if __name__ == '__main__':
    run()