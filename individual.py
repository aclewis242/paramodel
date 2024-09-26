from gen_funcs import *
from func_lib import *
from allele import *
import random
import time

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
    all_trans_bias: dict[str, float] = {}
    gtf_wgts: dict[str, float] = {}
    num_genes: int=0

    def __init__(self, gnts: list[str]=[], gnt: str='', gdm=wf, tps: list[list[float]]=[], gr: list[int]=[], rng: np.random.Generator=None,
                  **kwargs):
        '''
        Initialises the individual.

        ### Parameters
        - `gnts`: The list of genotypes that are allowed to be present within the individual
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
        self.num_genes = self.pc*self.ploidy
        self.genotype_freqs = dict.fromkeys(gnts, 0)
        if gnt: self.genotype_freqs[gnt] = self.pc
        if self.is_dip:
            if not gr: self.gene_range: list[int] = list(range(self.num_genes+1))
            else: self.gene_range = gr
            if not tps: self.trans_ps = gdm(self.num_genes)
            else: self.trans_ps = tps
        else:
            self.trans_ps = [[]]
            self.gene_range: list[int] = []
        if rng is None: self.rng = np.random.default_rng()
        else: self.rng = rng
        if self.pc_to_transmit > self.pc: self.pc_to_transmit = self.pc
        self.marked_for_death = False
        self.is_new_inf = True

    def simPara(self, times: list):
        '''
        Simulates the parasites' behavior. (`times` is for speed-recording purposes; refer to sim_lib for details.)
        '''
        for i in range(self.para_gens):
            # tm = time.time()
            if self.do_mutation: times = self.mutate(times)
            # self.checkGtfs('mutate')
            # times[13] += time.time() - tm
            times = self.genDrift(times)
            # self.checkGtfs('gendrift')
            if self.do_sr: self.genotype_freqs = self.reproduce()
        return times
    
    def genDrift(self, times: list):
        '''
        Simulates genetic drift & selection. (`times` is for speed-recording purposes; refer to sim_lib for details.)
        '''
        if self.file is not None: self.file.write('\t'.join([str(self.genotype_freqs[g]) for g in self.genotype_freqs])+'\n')
        if not self.is_mixed: return times
        tm = time.time()
        allele_freqs = self.getAlleleFreqs()
        times[6] += time.time() - tm
        tm = time.time()
        new_genotypes = []
        if self.num_alleles != 1: new_genotypes = ['']*self.pc
        curr_ind = len(allele_freqs)
        tba = ''
        # new_allele_freqs = allele_freqs.copy()
        times[7] += time.time() - tm
        for a in allele_freqs:
            curr_ind -= 1
            all_freq = allele_freqs[a]
            if not all_freq or all_freq == self.num_genes: continue
            tm = time.time()
            # if self.is_dip:
            #     print(f'omg something happening?! gtfs {self.genotype_freqs}, sel bias {self.all_sel_bias}')
            #     exit()
            all_prop = all_freq/self.num_genes
            asb = self.all_sel_bias[a]
            w_avg = asb*all_prop + (1 - all_prop)
            all_prop *= asb/w_avg
            # potential speed increase: condition on new_allele_freqs not being a nothing (0 or 1)
            # if self.is_dip and (new_allele_freqs[a] and new_allele_freqs[a] != self.num_genes):
            #     print(f'omg something happening?! gtfs {self.genotype_freqs}, sel bias {self.all_sel_bias}')
            #     print(new_allele_freqs)
            #     exit()
            times[8] += time.time() - tm
            tm = time.time()
            if self.is_dip:
                all_trans = round(all_prop*self.num_genes)
                if not all_trans or all_trans == self.num_genes: all_prop = float(bool(all_trans))
                else: all_prop = random.choices(self.gene_range, self.trans_ps[all_trans])[0]/self.num_genes
            j = 0
            times[15] += time.time() - tm
            if not all_prop or all_prop == 1:
                tm = time.time()
                all_prop_bool = bool(all_prop) # note: consider conditioning on hap instead of using ploidy (speed?)
                self.genotype_freqs[self.ploidy*a] = all_prop_bool*self.pc
                self.genotype_freqs[self.ploidy*chr(ord(a)+32)] = (not all_prop_bool)*self.pc
                if self.is_dip: self.genotype_freqs[a + chr(ord(a)+32)] = 0
                times[12] += time.time() - tm
                continue
            if self.is_hap and self.num_alleles == 1:
                tm = time.time()
                pc_flt = float(self.pc)
                gtfs_big = round(all_prop*pc_flt)
                gtfs_sml = round((1-all_prop)*pc_flt)
                if gtfs_big + gtfs_sml > pc_flt: gtfs_sml = self.pc - gtfs_big
                self.genotype_freqs[a] = gtfs_big
                self.genotype_freqs[chr(ord(a)+32)] = gtfs_sml
                times[9] += time.time() - tm
            else:
                tm = time.time()
                probs = self.ploidyProbs(all_prop)
                all_dist = self.rng.multinomial(n=self.pc, pvals=probs)
                times[10] += time.time() - tm
                if self.num_alleles == 1:
                    tm = time.time()
                    self.genotype_freqs[2*chr(ord(a)+32)] = all_dist[0]
                    self.genotype_freqs[a + chr(ord(a)+32)] = all_dist[1]
                    self.genotype_freqs[2*a] = all_dist[2]
                    times[11] += time.time() - tm
                else: # only used for multi-locus case, which (at least right now) should never arise. will greatly reduce speed
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
                    tm = time.time()
                    self.genotype_freqs = self.genotype_freqs.fromkeys(self.genotype_freqs, 0)
                    for n_g in new_genotypes: self.genotype_freqs[n_g] += 1
                    times[14] += time.time() - tm
        # tm = time.time()
        # self.storeData()
        # times[12] += time.time() - tm
        return times

    def mutate_old(self):
        '''
        Effects the mutations observed over a single generation. Based on an older and slower algorithm -- more general than the now-standard
        one, but much slower. Used only if there are mutations in a diploid, multi-locus individual (not currently part of the simulation)
        '''
        num_muts = self.rng.poisson(self.pc*self.num_alleles*self.mut_chance)
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
    
    def mutate(self, times: list):
        '''
        Effects the mutations observed over a single generation.
        '''
        if self.num_alleles != 1 or self.is_dip:
            self.mutate_old()
            return times
        tm = time.time()
        param_base = self.pc*self.num_alleles*self.mut_chance
        pre_gtfs = self.genotype_freqs.copy()
        times[1] += time.time() - tm
        for gt in self.genotype_freqs:
            tm = time.time()
            freq = pre_gtfs[gt]/self.pc
            times[4] += time.time() - tm
            if not freq: continue
            tm = time.time()
            num_muts = self.rng.poisson(param_base*freq)
            times[13] += time.time() - tm
            if not num_muts: continue
            tm = time.time()
            self.genotype_freqs[gt] -= num_muts
            self.genotype_freqs[gt.swapcase()] += num_muts
            times[14] += time.time() - tm
        return times
    
    def reproduce(self, s_d: dict[str, int]={}):
        '''
        Models sexual reproduction based on the given strain distribution. Note: haploid input distributions are okay, but it necessarily
        returns a diploid output distribution!

        (Currently deprecated)
        '''
        print('sr')
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

    def getAlleleFreqs(self):
        '''
        Gets the frequencies of each allele in the genotypes present in the individual's parasites. Used primarily in `genDrift`.
        '''
        rv: dict[str, int] = {}
        if self.num_alleles != 1:
            for g in self.genotype_freqs:
                genes = g.split('.')
                for gn in genes:
                    locus = gn[0].upper() # records mutated (uppercase) allele explicitly, wild (lowercase) implicitly
                    count = gn.count(locus)
                    if locus not in rv: rv[locus] = count*self.genotype_freqs[g]
                    else: rv[locus] += count*self.genotype_freqs[g]
        else:
            keys = list(self.genotype_freqs.keys())
            # Note: assumes genotypes are ordered from more uppercase to less uppercase (D-d, DD-Dd-dd)!
            keys_0 = keys[0]
            if self.is_hap: rv[keys_0] = self.genotype_freqs[keys_0]
            else: rv[keys_0[0]] = 2*self.genotype_freqs[keys_0] + self.genotype_freqs[keys[1]]
        return rv

    def getGenotypes(self):
        '''
        Returns all present genotypes as a list.
        '''
        return [gt for gt in self.genotype_freqs if self.genotype_freqs[gt]]
    
    def getGenotypeTransWeights(self):
        '''
        Get the genotypes' transmission weights.

        (Note: Intended to be used with transmission biases, but those are not currently implemented, so it just uses the base frequencies.)
        '''
        # return normalise(np.array(list(self.genotype_freqs.values())))
        return [self.genotype_freqs[gt]/self.pc for gt in self.genotype_freqs]

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
    
    def infectSelf(self, pc_num: int, strn: str, do_return: bool=False) -> list[str]:
        '''
        Infects the individual with the given strain using the given number of parasites.
        '''
        # if self.file is not None: self.file.write('\tinfectSelf\n')
        if pc_num > self.pc: pc_num = self.pc
        old_gnts = self.getGenotypes()
        to_replace = self.infectMult(pc_num)
        strn_match = self.match(strn)
        for stn in to_replace:
            self.genotype_freqs[stn] -= to_replace[stn]
            self.genotype_freqs[strn_match] += to_replace[stn]
        # self.checkGtfs('infectSelf')
        if do_return: return list(set(old_gnts) - set(self.getGenotypes()))
    
    def infectSelfMult(self, mix: dict[str, int]):
        '''
        Infects the individual with the given strain distribution.
        '''
        # if self.file is not None: self.file.write('\tinfectSelfMult\n')
        if self.do_sr:
            self.is_new_inf = True
            mix = self.reproduce(mix)
        mix_sum = sum(mix.values())
        if mix_sum >= self.pc:
            self.setToMix(mix, mix_sum=mix_sum)
            return
        for strn in mix: self.infectSelf(mix[strn], strn)
    
    def setToMix(self, mix: dict[str, int], mix_sum: int=0):
        '''
        Sets the individual's parasite distribution to the given strain distribution.
        '''
        # if self.file is not None: self.file.write('\tsetToMix\n')
        pc_transmitted = 0
        if not mix_sum: pc_transmitted = sum(mix.values())
        else: pc_transmitted = mix_sum
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
        # self.checkGtfs('s2mix')

    def match(self, s2m: str):
        '''
        Matches the given strain to the format required of this individual type (i.e., haploid or diploid).
        '''
        m2u = dipify
        if self.is_hap: m2u = hapify
        return m2u(s2m)
    
    def matchDist(self, sd2m: dict[str, int]):
        '''
        Takes a given strain distribution and matches it to the format required of this individual (i.e. hap/diploid).
        '''
        new_dist = dict.fromkeys(self.genotype_freqs, 0)
        for s in sd2m: new_dist[self.match(s)] += sd2m[s]
        return new_dist
    
    def correction_det(self, sn: str=''):
        '''
        'Corrects' for multi-counting by weighting the individual according to either the number of strains present within it (no parameter),
        or how prevalent the given strain is.
        '''
        if sn: return self.genotype_freqs[sn]/self.pc
        else: return 1/len(self.getGenotypes())
    
    def correction(self, sn: str=''):
        '''
        Returns either true or false with a probability of 1/(num. strains present) for no parameter, or the frequency of the given strain.
        '''
        return random.random() < self.correction_det(sn)
    
    def ploidyProbs(self, a_p: float=0.5):
        if self.is_hap: return [1-a_p, a_p]
        else: return [(1-a_p)**2, 2*a_p*(1-a_p), a_p**2]

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
    
    def checkGtfs(self, loc: str=''):
        if sum(self.genotype_freqs.values()) != self.pc:
            print(f'({loc}) pc {self.pc}, gtfs {self.genotype_freqs} (sum {sum(self.genotype_freqs.values())})')
            exit()
    
    @property
    def is_dip(self):
        return not self.is_hap
    
    @is_dip.setter
    def is_dip(self, value: bool):
        self.is_hap = not value
    
    @property
    def ploidy(self):
        return self.is_dip + 1
    
    # @property
    # def num_genes(self):
    #     return self.pc*self.ploidy
    
    @property
    def is_mixed(self):
        num_gnts = -1
        for gt in self.genotype_freqs:
            if self.genotype_freqs[gt]: num_gnts += 1
        return bool(num_gnts)
    
    def __str__(self):
        return f'{self.genotype_freqs}'
    
    def __repr__(self):
        return self.__str__()