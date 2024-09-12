from gen_funcs import *
from func_lib import *
from allele import *
import random
import time

import matplotlib.pyplot as plt # temp

class individual:
    '''
    The class for explicitly-modelled individuals.
    '''
    pc = 0
    genotype_freqs: dict[str, int] = {}
    trans_ps = [[]]
    file = None
    is_hap = False
    do_sr = False
    mut_chance = 0.0
    para_gens = 1
    do_mutation = False
    rng: np.random.Generator = None
    gene_range: list[int] = []
    pc_to_transmit = 0
    marked_for_death = False
    store_chance = 0.0
    is_new_inf = True
    num_alleles = -1
    all_sel_bias: dict[str, float] = {}
    bias_strength: 1.
    all_trans_bias: dict[str, float] = {}
    gtf_wgts: dict[str, float] = {}

    def __init__(self, alleles: list[allele]=[], gnt: str='', gdm=wf, tps: list[list[float]]=None, **kwargs):
        '''
        Initialises the individual.

        ### Parameters
        - `alleles`: The list of allele objects that the genotypes can be built from
        - `gnt`: The starting genotype of the parasites
        - `gdm`: The genetic drift model to use (ref. gen_funcs for options). Defaults to Wright-Fisher
        - `tps`: The matrix of transition probabilities to use for genetic drift. Generally the same as the output of `gdm`, but
            they can be pre-provided for the sake of speed
        - `pc`: The number of parasites present in the individual
        - `is_hap`: Whether or not the individual's parasites are haploid
        - `do_sr`: Whether or not sexual reproduction occurs inside the individual
        - `mut_chance`: The chance that a parasite will experience a mutation in an individual generation
        - `para_gens`: The number of generations the parasites go through per time step
        - `do_mutation`: Whether or not mutations can occur inside the individual
        - `pc_to_transmit`: The number of parasites to transmit during a mixed infection
        - `store_chance`: The chance that the individual's parasites' genetic makeup will be recorded in a file
        '''
        self.__dict__.update(kwargs)
        self.genotype_freqs: dict[str, int] = {}
        self.num_alleles = len(alleles)
        gnts = genGenotypes(alleles, self.is_hap)
        for g in gnts: self.genotype_freqs[g] = 0
        if gnt: self.genotype_freqs[gnt] = self.pc
        if tps is None: self.trans_ps = gdm(self.num_genes)
        else: self.trans_ps = tps
        self.rng = np.random.default_rng()
        self.gene_range: list[int] = list(range(self.num_genes+1))
        if self.pc_to_transmit > self.pc: self.pc_to_transmit = self.pc
        self.marked_for_death = False
        self.is_new_inf = True
        # self.all_trans_bias: dict[str, float] = {} # are the all_sel_biases at diff memory addresses?

    def simPara(self, times: list):
        '''
        Simulates the parasites' behavior. (`times` is for speed-recording purposes; refer to sim_lib for details.)
        '''
        for i in range(self.para_gens):
            if self.do_mutation: self.mutate()
            times = self.genDrift(times)
            if self.do_sr: self.genotype_freqs = self.reproduce()
        return times
    
    def genDrift(self, times: list):
        if self.file is not None: self.file.write('\t'.join([str(self.genotype_freqs[g]) for g in self.genotype_freqs])+'\n')
        tm = time.time()
        allele_freqs = self.getAlleleFreqs()
        times[9] += time.time() - tm # kind of expensive
        new_allele_freqs = allele_freqs.copy()
        new_genotypes = ['']*self.pc
        curr_ind = len(allele_freqs)
        tba = ''
        for a in allele_freqs:
            curr_ind -= 1
            tm = time.time()
            new_allele_freqs[a] = random.choices(self.gene_range, self.trans_ps[allele_freqs[a]])[0]
            times[6] += time.time() - tm
            tm = time.time()
            all_prop = self.bias(new_allele_freqs[a]/self.num_genes, self.all_sel_bias[a])
            # xs = np.array(range(11))/10
            # ys = [self.bias(x, a) for x in xs]
            # # print(self.genotype_freqs)
            # print(self.all_sel_bias)
            # plt.plot(xs, ys)
            # plt.show()
            # exit()
            probs = self.ploidyProbs(all_prop)
            times[7] += time.time() - tm
            tm = time.time()
            all_dist = self.rng.multinomial(n=self.pc, pvals=probs)
            times[10] += time.time() - tm # equally kind of expensive
            j = 0
            tm = time.time()
            for a_d_i in range(self.ploidy+1): # essentially a pc-length for loop
                for k in range(all_dist[a_d_i]):
                    if a_d_i == 0: tba = self.ploidy*chr(ord(a)+32) # changes char to lowercase quicker than .lower() method
                    elif a_d_i == 1:
                        if self.is_hap: tba = a
                        else: tba = (a + chr(ord(a)+32))
                    elif a_d_i == 2: tba = 2*a
                    if curr_ind: tba += '.'
                    new_genotypes[j] += tba
                    j += 1
            times[11] += time.time() - tm # most expensive but not by a huge amount
            tm = time.time()
            if curr_ind: random.shuffle(new_genotypes)
            times[12] += time.time() - tm # second most expensive
        self.genotype_freqs = self.genotype_freqs.fromkeys(self.genotype_freqs, 0)
        tm = time.time()
        for n_g in new_genotypes: self.genotype_freqs[n_g] += 1
        times[14] += time.time() - tm
        tm = time.time()
        times[13] += time.time() - tm
        self.storeData()
        return times

    def mutate(self):
        '''
        Effects the mutations observed over a single generation.
        '''
        num_muts = self.rng.binomial(self.pc*self.num_alleles, self.mut_chance)
        for i in range(num_muts):
            mut_src = self.infect()
            gnt_split = mut_src.split('.')
            all_idx = random.choice(list(range(self.num_alleles)))
            all_to_mut = gnt_split[all_idx]
            all_chances = dictify(genAlleles(all_to_mut), self.ploidyProbs())
            del all_chances[all_to_mut]
            new_sum = sum(all_chances.values())
            all_chances = {al: all_chances[al]/new_sum for al in all_chances}
            new_all = random.choices(list(all_chances.keys()), list(all_chances.values()))[0]
            gnt_split[all_idx] = new_all
            mut_tgt = '.'.join(gnt_split)
            self.genotype_freqs[mut_src] -= 1
            self.genotype_freqs[mut_tgt] += 1
            if self.file is not None: self.file.write('\tmut\n')
    
    def reproduce(self, s_d: dict[str, int]={}):
        '''
        Models sexual reproduction based on the given strain distribution. Note: haploid input distributions are okay, but it necessarily
        returns a diploid output distribution!
        '''
        if not s_d: s_d = self.genotype_freqs
        if (not self.is_new_inf or self.is_hap) or not self.do_sr: return s_d
        new_dist = self.matchDist(dict.fromkeys(s_d, 0))
        total_num = sum(s_d.values())
        for i in range(total_num):
            parents = [self.infect(s_d).split('.') for j in range(2)]
            new_gnt = ''
            for j in range(self.num_alleles):
                new_all = ''
                for par in parents: new_all += random.choice(par[j])
                new_all_lst = list(new_all)
                new_all_lst.sort()
                new_gnt += (''.join(new_all_lst) + '.')
            new_gnt = new_gnt[:-1]
            new_dist[new_gnt] += 1
        self.is_new_inf = False
        if self.file is not None:
            self.file.write(f'orig: {self.genotype_freqs}\n')
            self.file.write(f'sr: {s_d} -> {new_dist}\n')
        return new_dist

    def getGenotypes(self):
        '''
        Returns all present genotypes as a list.
        '''
        return [gt for gt in self.genotype_freqs if self.genotype_freqs[gt]]
    
    def getGenotypeTransWeights(self) -> np.ndarray[float]:
        # print(self.gtf_wgts)
        # print(self.genotype_freqs)
        def normDictVals(dct: dict):
            return normalise(np.array(list(dct.values())))
        gtf_vals = normDictVals(self.genotype_freqs)
        gtfs_norm = dictify(self.genotype_freqs.keys(), gtf_vals)
        gtfs_wgt = {gt: gtfs_norm[gt]**(self.gtf_wgts[gt]**self.bias_strength) for gt in gtfs_norm}
        # exit()
        return normDictVals(gtfs_wgt)

    def infectMult(self, num: int=1) -> dict[str, int]:
        '''
        Performs multiple infections. Returns a dict of strain to # times infected.
        '''
        return dictify(self.genotype_freqs.keys(), self.rng.multinomial(num, self.getGenotypeTransWeights()))

    def infect(self, gtfs: dict[str, int]={}):
        '''
        Performs an infection.
        '''
        if not gtfs: gtfs = self.genotype_freqs
        return random.choices(list(gtfs.keys()), self.getGenotypeTransWeights())[0]
    
    def getAlleleFreqs(self):
        '''
        Gets the frequencies of each allele in the genotypes present in the individual's parasites. Used primarily for genetic drift.
        '''
        rv: dict[str, int] = {}
        for g in self.genotype_freqs:
            genes = g.split('.')
            for gn in genes:
                locus = gn[0].upper() # records 'dominant' (or just capital allele in haploid case) explicitly, 'recessive' implicitly
                count = gn.count(locus)
                if locus not in rv: rv[locus] = count*self.genotype_freqs[g]
                else: rv[locus] += count*self.genotype_freqs[g]
        return rv
    
    def infectSelf(self, pc_num: int, strn: str):
        '''
        Infects the individual with the given strain using the given number of parasites.
        '''
        if self.file is not None: self.file.write('\tinfectSelf\n')
        if pc_num > self.pc: pc_num = self.pc
        old_gnts = self.getGenotypes()
        to_replace = self.infectMult(pc_num)
        for stn in to_replace:
            self.genotype_freqs[stn] -= to_replace[stn]
            self.genotype_freqs[self.match(strn)] += to_replace[stn]
        if sum(self.genotype_freqs.values()) != self.pc:
            print(f'(inside indv) pc {self.pc}, gtfs {self.genotype_freqs}. submitted strain {strn}, pc_num {pc_num}')
            exit()
        return list(set(old_gnts) - set(self.getGenotypes()))
    
    def infectSelfMult(self, mix: dict[str, int]):
        '''
        Infects the individual with the given strain distribution.
        '''
        if self.file is not None: self.file.write('\tinfectSelfMult\n')
        if self.do_sr:
            self.is_new_inf = True
            mix = self.reproduce(mix)
        if sum(mix.values()) >= self.pc:
            self.setToMix(mix)
            return
        for strn in mix: self.infectSelf(mix[strn], strn)
    
    def setToMix(self, mix: dict[str, int]):
        '''
        Sets the individual's parasite distribution to the given strain distribution.
        '''
        if self.file is not None: self.file.write('\tsetToMix\n')
        pc_transmitted = sum(mix.values())
        rem = 0.
        matched = ''
        self.genotype_freqs = dict.fromkeys(self.genotype_freqs, 0)
        for strn in mix:
            matched = self.match(strn)
            amt_raw = self.pc*(mix[strn]/pc_transmitted) + rem
            amt = int(amt_raw)
            rem = amt_raw - amt
            self.genotype_freqs[matched] += amt
        if rem > 0.999: self.genotype_freqs[matched] += 1
        if sum(self.genotype_freqs.values()) != self.pc:
            print(f'(inside s2mix) pc {self.pc}, gtfs {self.genotype_freqs}. submitted mix {mix}')
            exit()

    def match(self, s2m: str):
        '''
        Matches the given strain to the format required of this individual type (i.e., haploid or diploid).
        '''
        m2u = dipify
        if self.is_hap: m2u = hapify
        return m2u(s2m)
    
    def matchDist(self, sd2m: dict[str, int]):
        new_dist = dict.fromkeys(self.genotype_freqs, 0)
        for s in sd2m: new_dist[self.match(s)] += sd2m[s]
        return new_dist
    
    def correction(self, sn: str=''):
        '''
        'Corrects' for overcounting by rolling a die, with success weighted either by the given strain's frequency or (if blank)
        the number of strains present.
        '''
        return True
        if sn: return random.random() < self.genotype_freqs[sn]/self.pc
        return random.random() < 1/len(self.getGenotypes())
        # num = 0
        # if sn: num = self.genotype_freqs[sn]/self.pc
        # else: num = 1/len(self.getGenotypes())
        # f = open('anti_stoch.dat')
        # num += float(f.readline())
        # f.close()
        # f = open('anti_stoch.dat', 'w')
        # f.write(str(num%1))
        # return num >= 1
    
    def correction_det(self, sn: str=''):
        # return 1
        if sn: return self.genotype_freqs[sn]/self.pc
        else: return 1/len(self.getGenotypes())
    
    def ploidyProbs(self, a_p: float=0.5):
        if self.is_hap: return [1-a_p, a_p]
        else: return [(1-a_p)**2, 2*a_p*(1-a_p), a_p**2]

    def bias(self, a_p: float, bias_num: float):
        return a_p**((2*(1-bias_num))**self.bias_strength)

    def storeData(self, force: bool=False):
        '''
        Stores genotype frequency data, maybe (depends on `store_chance`).
        '''
        if self.file is None and (random.random() <= self.store_chance or force):
            self.file = open(f'{int(random.random()*1e6)}.dat', 'x')
            [self.file.write(f'{gnt}\t') for gnt in self.genotype_freqs]
            self.file.write('\n')
    
    def rebuild(self):
        '''
        Currently deprecated; may eventually be used for deep copying.
        '''
        return
    
    @property
    def is_dip(self):
        return not self.is_hap
    
    @is_dip.setter
    def is_dip(self, value: bool):
        self.is_hap = not value
    
    @property
    def ploidy(self):
        return self.is_dip + 1
    
    @property
    def num_genes(self):
        return self.pc*self.ploidy
    
    def __str__(self):
        hapdip = 'di'
        if self.is_hap: hapdip = 'ha'
        return f'{hapdip}ploid'
    
    def __repr__(self):
        return self.__str__()