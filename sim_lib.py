from model import *
from allele import *
import numpy as np
import time

def simShell(tmax: float, mdls: list[SIR], nt: float=2e5, alleles: list[allele]=[]):
    '''
    Manages the time iterations of the simulation.

    ### Parameters
    - `tmax`: The maximum amount of time to run the simulation for.
    - `mdls`: A list containing the models governing each population (initialised with parameters & populations).
    - `nt`: The number of time steps to use, as a 'float' (e.g. 2e5 - integer in floating-point form).
    - `alleles`: The list of all possible alleles (allele objects). Irrelevant if not using the allele model.
    - `is_haploid`: Whether the model is haploid or diploid (bool).

    ### Returns
    - `ts`: A NumPy array of the times visited by the simulation (not equispaced nor integer).
    - `ps`: A NumPy array that contains the populations (flattened) at each time. Rows index population, columns index time.
    - `times`: An array containing the computation times of the various components of the method.
    - `pops`: A list of all population objects.
    '''
    dt = tmax/(nt - 1)
    nt = int(nt)

    vec_strain_2_mdl = {}
    is_haploid = False
    for m in mdls:
        if m.pop.is_hap: is_haploid = True
    is_hyb = False
    for m in mdls:
        if m.pop.is_dip and is_haploid: is_hyb = True
    if len(alleles):
        genotypes = genGenotypes(alleles, is_haploid)
        if is_hyb: genotypes_dip = genGenotypes(alleles, False)
        else: genotypes_dip = []
        new_models = []
        mdls.sort(key=lambda x: x.is_vector, reverse=True)
        for m in mdls: # vec mdl is assumed to be first
            gt = genotypes
            if m.is_vector and is_hyb: gt = genotypes_dip
            gt.reverse() # to put 'recessive' (small) first, as it's assumed to be the 'wild' type
            for g in gt:
                vec_mdls = [None]
                new_models_temp = []
                if not m.is_vector: vec_mdls = vec_strain_2_mdl[g]
                new_models_temp = [m.updateGenotype(g, alleles, vec_mdl) for vec_mdl in vec_mdls]
                new_models_temp = new_models_temp[:1]
                for new_model in new_models_temp:
                    if 'init' in new_model.pop.inf.keys():
                        new_model.pop.addPop(new_model.pop.getPop(), g)
                        del new_model.pop.inf['init']
                        del new_model.pop.rec['init']
                    new_models += [new_model]
                    if m.is_vector:
                        g2 = g
                        if is_hyb: g2 = hapify(g)
                        if g2 in vec_strain_2_mdl.keys(): vec_strain_2_mdl[g2] += [new_model]
                        else: vec_strain_2_mdl[g2] = [new_model]
        mdls = new_models

    [m.setRs() for m in mdls]
    strain_2_mdl: dict[str, list[SIR]] = {}
    for m in mdls:
        if m.sn not in strain_2_mdl.keys(): strain_2_mdl[m.sn] = [m]
        else: strain_2_mdl[m.sn] += [m]
    for s in strain_2_mdl:
        s_mdls = strain_2_mdl[s]
        vec_mdl = [sm for sm in s_mdls if sm.is_vector]
        if len(vec_mdl):
            vec_mdl = vec_mdl[0]
            if not is_hyb: s_mdls.pop(s_mdls.index(vec_mdl))
            else: continue
        else: vec_mdl = strain_2_mdl[dipify(s)][0]
        s_mdls_2: list[SIR] = []
        for sm in s_mdls:
            if sm not in s_mdls_2: s_mdls_2 += [sm]
        s_mdls = s_mdls_2
        [print(f'{sm} r0: {sm.r0(vec_mdl)}') for sm in s_mdls]
    ts_i = np.array(range(int(nt)))
    ps_init = np.empty(shape=(nt, len(mdls), len(mdls[0].pop.getAllPop())))
    ps = listify(ps_init)
    pops = [m.pop for m in mdls]
    pops_check = []
    for p in list(pops):
        if hex(id(p)) not in pops_check: pops_check += [hex(id(p))]
        else: pops.remove(p)
    num_pops = len(pops)
    times = [0 for i in range(16)]
    num_Rs = len(mdls[0].Rs)
    num_mdls = len(mdls)
    all_Rs = np.array([0.0 for i in range(num_mdls*num_Rs)])
    vec_pop = [p for p in pops if p.is_vector][0]
    for p in pops:
        max_inf = -1
        max_strn = ''
        for s in p.inf:
            if p.inf[s] > max_inf:
                max_inf = p.inf[s]
                max_strn = s
        p.inf[max_strn] = 0
        p.inf[max_strn.lower()] = max_inf

    # [print(m.__dict__) for m in mdls]
    # exit()

    for i in ts_i:
        tm = time.time()
        for j in range(num_mdls):
            mdls[j].setRs()
            for k in range(num_Rs):
                all_Rs[j*num_Rs+k] = mdls[j].Rs[k]
        times[0] += time.time() - tm # a little expensive. scales poorly with simulation complexity
        tm = time.time()
        sum_Rs = sum(all_Rs)
        times[1] += time.time() - tm
        tm = time.time()
        if not sum_Rs: break
        times[2] += time.time() - tm
        tm = time.time()
        Xs = adaptSim(all_Rs/sum_Rs, sum_Rs, dt)
        times[3] += time.time() - tm # most expensive
        for i_m in range(num_mdls):
            for i_r in range(num_Rs):
                tm = time.time()
                rpt = Xs[i_m*num_Rs+i_r]
                if rpt: mdls[i_m].trans(i_r, rpt)
                times[4] += time.time() - tm # 2nd most expensive
        tm = time.time()
        for i_p in range(num_pops): ps[i][i_p] = pops[i_p].getAllPop()
        times[5] += time.time() - tm
        # All remaining measurement blocks are left over from a more complicated and inefficient time
        # Preserved in case they prove useful sometime in the future

        tm = time.time()
        times[6] += time.time() - tm
        tm = time.time()
        for p in pops:
            alive: list[individual] = []
            for indv_i in range(len(p.individuals)):
                indv = p.individuals[indv_i]
                if not indv.marked_for_death:
                    init_gts = indv.getGenotypes()
                    times = indv.simPara(times)
                    fin_gts = indv.getGenotypes()
                    if fin_gts != init_gts:
                        gts_rmv = list(set(init_gts) - set(fin_gts))
                        gts_add = list(set(fin_gts) - set(init_gts))
                        for gt in gts_rmv: p.inf[gt] -= 1
                        for gt in gts_add: p.inf[gt] += 1
                    alive += [indv]
            p.individuals = alive
        times[15] += time.time() - tm
        tm = time.time()
        print(f'{int(100*i/nt)}%; vec indvs: {len(vec_pop.individuals)} ', end='\r')
        times[8] += time.time() - tm
        tm = time.time()
        times[9] += time.time() - tm
        tm = time.time()
        times[10] += time.time() - tm
        tm = time.time()
        times[11] += time.time() - tm
        tm = time.time()
        # [print(p.__dict__) for p in pops]
        times[12] += time.time() - tm
        tm = time.time()
        times[13] += time.time() - tm
        tm = time.time()
        times[14] += time.time() - tm
    times[15] -= (sum(times[9:15]) + sum(times[6:8]))
    return ts_i*dt, ps, times, pops

def adaptSim(ps: np.ndarray[float], sum_Rs: float, dt: float):
    '''
    Adaptively picks the best way to estimate the results of the model. Returns an array containing the number of times each event
    occurs (ordered the same way as the given probabilities).

    ### Parameters
    - `ps`: The relative probabilities of each event, as a NumPy array.
    - `sum_Rs`: The net rate of all events.
    - `dt`: The size of the time step.
    '''
    Xs = 0*ps
    p_cond = 0
    rng = np.random.default_rng()
    N = int(sum_Rs*dt) # to save on time, the random variable this should technically be has been replaced with its avg value
    det_thres = 0.1
    for i in range(len(ps)):
        if p_cond >= 1 or N <= 0: break
        p = ps[i]/(1 - p_cond)
        if p > 1: p = 1
        if p <= 0: continue
        if False: # re-enable this if speed becomes particularly necessary again (it makes ensuring accuracy... annoying at best)
            if N > 200 and (p > det_thres and p < 1-det_thres): Xs[i] = int(N*p)              # Deterministic case
            elif N > 1000:
                if N*p < 25 or N*(1-p) < 25: Xs[i] = rng.poisson(lam=N*p)       # Large-ish N, p close to 0 or 1
                else: Xs[i] = abs(int(rng.normal(loc=N*p, scale=N*p*(1-p))))    # Large-ish N, p close to neither 0 nor 1
            else: Xs[i] = rng.binomial(n=N, p=p)                                # Small N
        else: Xs[i] = rng.binomial(n=N, p=p)
        N -= Xs[i]
        p_cond += ps[i]
    return Xs