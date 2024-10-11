from model import *
from allele import *
from time import time
import numpy as np

def simShell(tmax: float, mdls: list[SIR], nt: float=2e5, alleles: list[allele]=[], weight_infs: bool=False, do_mix_start: bool=False,
             num_hist: int=5):
    '''
    Manages the time iterations of the simulation.

    ### Parameters
    - `tmax`: The maximum amount of time to run the simulation for.
    - `mdls`: A list containing the models governing each population (initialised with parameters & populations).
    - `nt`: The number of time steps to use, as a 'float' (e.g. 2e5 - integer in floating-point form).
    - `alleles`: The list of all possible alleles (allele objects). Irrelevant if not using the allele model.
    - `weight_infs`: Whether or not to weight infected data according to individuals' genotype frequencies.
    - `do_mix_start`: Whether or not to have a mixed distribution of infected individuals (wild & mutated) or uniform (just wild).
    - `num_hist`: The number of histograms to generate (varying over time).

    ### Returns
    - `ts`: A NumPy array of the times visited by the simulation (neither equispaced nor integer).
    - `ps`: A NumPy array that contains the populations (flattened) at each time. Rows index population, columns index time.
    - `times`: A list containing the computation times of the various components of the method.
    - `pops`: A list of all population objects.
    - `ps_unwgt`: Like `ps`, but with unweighted infected information.
    - `vpis`: A list of the total number of infected vectors at each time.
    - `hpis`: Like `vpis`, but for hosts instead of vectors.
    - `hists_v`: Histogram data for vectors.
    - `hists_h`: Histogram data for hosts.
    - `hist_tms`: The times at which histogram data were taken.
    '''
    # tm = time()
    times = [0 for i in range(21)]
    dt = tmax/(nt - 1)
    nt = int(nt)
    if num_hist < 2 and num_hist != 0: num_hist = 2
    hist_fac = ((num_hist - 1)/nt)*bool(num_hist)

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
            m.pop.gnts = gt.copy()
            gt.reverse() # to put 'recessive' (small) first, as it's assumed to be the 'wild' type
            for g in gt:
                vec_mdls = [None]
                new_models_temp = []
                if not m.is_vector: vec_mdls = vec_strain_2_mdl[g]
                new_models_temp = [m.updateGenotype(g, alleles, vec_mdl) for vec_mdl in vec_mdls]
                new_models_temp = new_models_temp[:1]
                for new_model in new_models_temp:
                    if 'init' in new_model.pop.inf.keys():
                        new_model.pop.addPop([0] + new_model.pop.getPop()[1:], g)
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
        # [print(f'{sm} r0: {sm.r0(vec_mdl)}') for sm in s_mdls]
    ts_i = np.array(range(int(nt)))
    ps_init = np.empty(shape=(nt, len(mdls), len(mdls[0].pop.getAllPop())))
    ps = listify(ps_init)
    ps_unwgt = listify(ps)
    pops = [m.pop for m in mdls]
    pops_check = []
    for p in list(pops):
        if hex(id(p)) not in pops_check: pops_check += [hex(id(p))]
        else: pops.remove(p)
    num_pops = len(pops)
    num_Rs = len(mdls[0].Rs)
    num_mdls = len(mdls)
    all_Rs = np.array([0.0 for i in range(num_mdls*num_Rs)])
    vec_pop = [p for p in pops if p.is_vector][0]
    host_pop = [p for p in pops if not p.is_vector][0]
    for p in pops:
        max_inf = -1
        max_strn = ''
        for s in p.inf:
            if p.inf[s] > max_inf:
                max_inf = p.inf[s]
                max_strn = s
        if do_mix_start:
            p.inf[max_strn] = max_inf/2
            p.inf[max_strn.lower()] = max_inf/2
            for ind in p.individuals[:int(len(p.individuals)/2)]:
                ind.genotype_freqs[max_strn.upper()] = ind.pc
                ind.genotype_freqs[max_strn.lower()] = 0
            shuffle(p.individuals)
            # [print(ind.genotype_freqs) for ind in p.individuals]
        else:
            p.inf[max_strn] = 0
            p.inf[max_strn.lower()] = max_inf
    [p.updateSelBiases(alleles) for p in pops]
    for p in pops:
        p.init_pop = p.tot_pop
        p.num_gnts = len(p.gnts)
    for m in mdls: m.bds = m.bd/m.pop.num_gnts

    hpis = []
    vpis = []
    hists_h = []
    hists_v = []
    hist_tms = []
    # times[0] += time() - tm
    for i in ts_i:
        # tm = time()
        for j in range(num_mdls):
            mdls[j].setRs()
            for k in range(num_Rs):
                all_Rs[j*num_Rs+k] = mdls[j].Rs[k]
        sum_Rs = sum(all_Rs)
        # times[1] += time() - tm
        # tm = time()
        Xs = adaptSim(all_Rs/sum_Rs, sum_Rs, dt)
        # times[3] += time() - tm
        # tm = time()
        for i_m in range(num_mdls):
            for i_r in range(num_Rs):
                rpt = Xs[i_m*num_Rs+i_r]
                if rpt: mdls[i_m].trans(i_r, rpt)
        # times[2] += time() - tm
        # tm = time()
        for i_p in range(num_pops):
            ps[i][i_p] = pops[i_p].getAllPop(weight=weight_infs)
            ps_unwgt[i][i_p] = pops[i_p].getAllPop()
        # times[3] += time() - tm
        for p in pops:
            # tm = time()
            # alive: list[individual] = []
            # shuffle(p.individuals)
            # times[4] += time() - tm
            for indv in p.individuals: indv.simPara()
                # # tm = time()
                # if sum(indv.genotype_freqs.values()) != indv.pc:
                #     print(f'indv with {indv.pc} parasites has gtfs {indv.genotype_freqs} (sum {sum(indv.genotype_freqs.values())})')
                #     exit()
                # # times[4] += time() - tm
                # if not indv.marked_for_death:
                    # tm = time()
                    # init_gts = indv.getGenotypes()
                    # times[5] += time() - tm
                    # indv.simPara()
                    # tm = time()
                    # fin_gts = indv.getGenotypes()
                    # times[5] += time() - tm
                    # tm = time()
                    # if fin_gts != init_gts:
                    #     # print(f'thing actually happening. init gts {init_gts}, fin gts {fin_gts}')
                    #     gts_rmv = list(set(init_gts) - set(fin_gts))
                    #     gts_add = list(set(fin_gts) - set(init_gts))
                    #     for gt in gts_rmv: p.inf[gt] -= 1
                    #     for gt in gts_add: p.inf[gt] += 1
                    #     for strn in p.inf:
                    #         if p.inf[strn] < 0: print('concern')
                    # alive += [indv]
                    # times[18] += time() - tm
            p.update()
            # p.individuals = alive
        # tm = time()
        vpi = len(vec_pop.individuals)
        hpi = len(host_pop.individuals)
        vpis += [vpi]
        hpis += [hpi]
        print(f'{int(100*i/nt)}%; vec indvs: {vpi}; host indvs: {hpi}; vec pop: {vec_pop.tot_pop}; host pop: {host_pop.tot_pop} ',
               end='\r')
        hist_check = hist_fac*i
        hist_check_2 = int(hist_fac*(i+1))
        hist_check_int = int(hist_check)
        if (hist_check_int == hist_check or hist_check_2 > hist_check_int) and num_hist:
            hists_h += [host_pop.getGntDist()]
            hists_v += [vec_pop.getGntDist()]
            hist_tms += [i*dt]
        # times[19] += time() - tm
    # tm = time()
    hists_h += [host_pop.getGntDist()]
    hists_v += [vec_pop.getGntDist()]
    hist_tms += [ts_i[-1]*dt]
    # times[20] += time() - tm
    return ts_i*dt, ps, times, pops, ps_unwgt, vpis, hpis, hists_v, hists_h, hist_tms

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