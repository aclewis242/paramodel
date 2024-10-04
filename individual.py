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
    is_hap = False
    mut_chance = 0.0
    para_gens = 1
    do_mutation = False
    rng: np.random.Generator = None
    gene_range: list[int] = []
    pc_to_transmit = 0
    marked_for_death = False
    all_sel_bias: dict[str, float] = {}
    all_transm_probs: dict[str, float] = {}
    num_genes: int = 0
    pc_flt: float = 0.0

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
        - `gr`: The range of allowed outcomes from the transition probability matrix (0-ploidy*pc). Not meant to vary, but taking it
            as an argument helps with speed
        - `rng`: The NumPy random number generator object to use. Again, not meant to vary, but making a new one for every individual
            would be unnecessary & time-consuming
        - `pc`: The number of parasites present in the individual
        - `is_hap`: Whether or not the individual's parasites are haploid
        - `mut_chance`: The chance that a parasite will experience a mutation in an individual generation
        - `para_gens`: The number of generations the parasites go through per time step
        - `do_mutation`: Whether or not mutations can occur inside the individual
        - `pc_to_transmit`: The number of parasites to transmit during a mixed infection
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

    def simPara(self, times: list):
        '''
        Simulates the parasites' behavior. (`times` is for speed-recording purposes; refer to sim_lib for details.)
        '''
        for i in range(self.para_gens):
            if self.do_mutation: times = self.mutate(times)
            times = self.genDrift(times)
        return times
    
    def genDrift(self, times: list):
        '''
        Simulates genetic drift & selection. (`times` is for speed-recording purposes; refer to sim_lib for details.)
        '''
        if not self.is_mixed: return times
        tm = time.time()
        allele_freqs = self.getAlleleFreqs()
        times[6] += time.time() - tm
        for a in allele_freqs:
            tm = time.time()
            all_freq = allele_freqs[a]
            times[7] += time.time() - tm
            if not all_freq or all_freq == self.num_genes: continue
            tm = time.time()
            all_prop = all_freq/self.num_genes
            asb = self.all_sel_bias[a]
            w_avg = asb*all_prop + (1. - all_prop)
            all_prop *= asb/w_avg
            times[8] += time.time() - tm
            tm = time.time()
            if self.is_dip:
                all_trans = round(all_prop*self.num_genes)
                # if not all_trans or all_trans == self.num_genes:
                #     all_prop = float(bool(all_trans))
                #     print(f'init cond: {all_prop}')
                # else:
                all_prop = random.choices(self.gene_range, self.trans_ps[all_trans])[0]/self.num_genes
            times[9] += time.time() - tm
            if not all_prop or all_prop == 1:
                tm = time.time()
                all_prop_bool = bool(all_prop) # note: consider conditioning on hap instead of using ploidy (speed?)
                self.genotype_freqs[self.ploidy*a] = all_prop_bool*self.pc
                self.genotype_freqs[self.ploidy*chr(ord(a)+32)] = (not all_prop_bool)*self.pc
                if self.is_dip: self.genotype_freqs[a + chr(ord(a)+32)] = 0 # switches to lowercase faster than .lower() method
                times[10] += time.time() - tm
                continue
            if self.is_hap:
                tm = time.time()
                gtfs_big = round(all_prop*self.pc_flt)
                gtfs_sml = round((1-all_prop)*self.pc_flt)
                if gtfs_big + gtfs_sml > self.pc: gtfs_sml = self.pc - gtfs_big
                self.genotype_freqs[a] = gtfs_big
                self.genotype_freqs[chr(ord(a)+32)] = gtfs_sml
                times[11] += time.time() - tm
            else:
                tm = time.time()
                probs = [(1-all_prop)**2, 2*all_prop*(1-all_prop), all_prop**2]
                all_dist = self.rng.multinomial(n=self.pc, pvals=probs)
                times[12] += time.time() - tm
                tm = time.time()
                self.genotype_freqs[2*chr(ord(a)+32)] = all_dist[0]
                self.genotype_freqs[a + chr(ord(a)+32)] = all_dist[1]
                self.genotype_freqs[2*a] = all_dist[2]
                times[13] += time.time() - tm
        return times

    def mutate(self, times: list):
        '''
        Effects the mutations observed over a single generation.
        '''
        if self.is_dip: return times
        tm = time.time()
        param_base = self.pc_flt*self.mut_chance
        pre_gtfs = self.genotype_freqs.copy()
        times[14] += time.time() - tm
        for gt in self.genotype_freqs:
            tm = time.time()
            mut_param = param_base*pre_gtfs[gt]/self.pc_flt
            times[15] += time.time() - tm
            if mut_param < 0.05: continue
            tm = time.time()
            num_muts = self.rng.poisson(mut_param)
            # print(f'num_muts {num_muts}, param {mut_param}')
            times[16] += time.time() - tm
            if not num_muts: continue
            tm = time.time()
            self.genotype_freqs[gt] -= num_muts
            self.genotype_freqs[gt.swapcase()] += num_muts
            times[17] += time.time() - tm
        return times

    def getAlleleFreqs(self):
        '''
        Gets the frequencies of each allele in the genotypes present in the individual's parasites. Used primarily in `genDrift`.
        '''
        rv: dict[str, int] = {}
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
    
    # def getGenotypeTransWeights(self):
    #     '''
    #     Get the genotypes' transmission weights.
    #     '''
    #     weights_unnorm = [self.all_transm_probs[gt[0]]*self.genotype_freqs[gt]/self.pc_flt for gt in self.genotype_freqs]
    #     wgts_sum = sum(weights_unnorm)
    #     return [wgt_unnorm/wgts_sum for wgt_unnorm in weights_unnorm], wgts_sum
    
    def getGenotypeTransWeights_unwgt(self):
        return [self.genotype_freqs[gt]/self.pc_flt for gt in self.genotype_freqs]

    # def infectMult(self, num: int=1) -> dict[str, int]:
    #     '''
    #     Performs multiple infections. Returns a dict of strain to # times infected.
    #     '''
    #     gtf_wgts, num_adj = self.getGenotypeTransWeights()
    #     print(f'gtf_wgts: {gtf_wgts}, num_adj: {num_adj}, num: {num}, gtfs: {self.genotype_freqs}')
    #     exit()
    #     return dictify(self.genotype_freqs.keys(), self.rng.multinomial(round(num*num_adj), gtf_wgts))

    def doesContactTransmit(self):
        gtf_wgt_sum = 0
        for gt in self.genotype_freqs: gtf_wgt_sum += self.all_transm_probs[gt[0]]*self.genotype_freqs[gt]/self.pc_flt
        return random.random() < gtf_wgt_sum

    def infectMix(self, pc_num: int=1, do_test_contact: bool=True):
        if do_test_contact:
            if not self.doesContactTransmit(): return dict.fromkeys(self.genotype_freqs, 0)
        return dictify(self.genotype_freqs.keys(), self.rng.multinomial(pc_num, self.getGenotypeTransWeights_unwgt()))

    def infect(self, gtfs: dict[str, int]={}):
        '''
        Performs an infection. (Obsolete?)
        '''
        if not gtfs: gtfs = self.genotype_freqs
        print('doing infect method (this shouldn\'t be happening anywhere, I think)')
        return random.choices(list(gtfs.keys()), self.getGenotypeTransWeights_unwgt())[0]
    
    def infectSelf(self, pc_num: int, strn: str, do_return: bool=False) -> list[str]:
        '''
        Infects the individual itself.

        ### Parameters
        - `pc_num`: The number of parasites in the transmission.
        - `strn`: The strain to infect it with.
        - `do_return`: Whether or not to return the genotypes no longer present afterwards.
        '''
        if pc_num > self.pc: pc_num = self.pc # simplify w/assumptions about relative pc sizes?
        old_gnts = self.getGenotypes()
        to_replace = self.infectMix(pc_num, do_test_contact=False)
        strn_match = self.match(strn)
        for stn in to_replace:
            self.genotype_freqs[stn] -= to_replace[stn]
            self.genotype_freqs[strn_match] += to_replace[stn]
        if do_return: return list(set(old_gnts) - set(self.getGenotypes()))
    
    def infectSelfMult(self, mix: dict[str, int]):
        '''
        Infects the individual with the given strain distribution.
        '''
        mix_sum = sum(mix.values())
        if mix_sum >= self.pc:
            self.setToMix(mix, mix_sum=mix_sum)
            return
        for strn in mix: self.infectSelf(mix[strn], strn)
    
    def setToMix(self, mix: dict[str, int], mix_sum: int=0):
        '''
        Sets the individual's parasite distribution to the given strain distribution. (If the sum of the distribution is already known,
        it can be provided as a keyword argument.)
        '''
        pc_transmitted = 0
        if not mix_sum: pc_transmitted = sum(mix.values())
        else: pc_transmitted = mix_sum
        rem = 0.
        matched = ''
        self.genotype_freqs = dict.fromkeys(self.genotype_freqs, 0)
        for strn in mix:
            matched = self.match(strn)
            amt_raw = self.pc_flt*(mix[strn]/pc_transmitted) + rem
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
        if sn: return self.genotype_freqs[sn]/self.pc_flt
        else: return 1/len(self.getGenotypes())
    
    def correction(self, sn: str=''):
        '''
        Returns either true or false with a probability of 1/(num. strains present) for no parameter, or the frequency of the given strain.
        '''
        return random.random() < self.correction_det(sn)
    
    def checkGtfs(self, loc: str=''):
        '''
        Makes sure the genotype freq. values align with the total parasite #. Takes some identifying string (e.g. the source of the call)
        as an argument.
        '''
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