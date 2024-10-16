from gen_funcs import *
from func_lib import *
from random import choices, random

class individual:
    '''
    The class for explicitly-modelled individuals.
    '''
    pc = 0                                  # Parasite population within the individual
    genotype_freqs: dict[str, int] = {}     # Dict of genotypes to 'frequencies' (absolute #s)
    is_hap = False                          # Whether or not this individual's parasites are haploid
    mut_chance = 0.0                        # Chance of mutation per parasite per generation
    para_gens = 1                           # Number of parasite generations per time step
    do_mutation = False                     # Whether or not this individual's parasites are allowed to mutate
    rng: np.random.Generator = None         # NumPy random number generator object
    pc_to_transmit = 0                      # The number of parasites that get transmitted during a successful infection
    marked_for_death = False                # Marked true to indicate that it has died/recovered and needs to be removed
    all_sel_bias: dict[str, float] = {}     # Allele selection biases (for genetic drift)
    all_transm_probs: dict[str, float] = {} # Allele transmission probabilities
    num_genes: int = 0                      # The total number of alleles throughout the individual's parasite population
    pc_flt: float = 0.0                     # Parasite population as a float
    main_all_char: str = ''                 # Mutated (uppercase) allele character
    scnd_all_char: str = ''                 # Wild (lowercase) allele character
    is_mixed_vec: bool = False              # Whether or not this individual started out as a mixed-strain vector

    def __init__(self, gnts: list[str]=[], gnt: str='', rng: np.random.Generator=None, **kwargs):
        '''
        Initialises the individual.

        ### Parameters
        - `gnts`: The list of genotypes that are allowed to be present within the individual
        - `gnt`: The starting genotype of the parasites
        - `rng`: The NumPy random number generator object to use. Not meant to vary, but making a new one for every individual
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
        if rng is None: self.rng = np.random.default_rng()
        else: self.rng = rng
        if self.pc_to_transmit > self.pc: self.pc_to_transmit = self.pc

    def simPara(self):
        '''
        Simulates the parasites' behavior.
        '''
        for i in range(self.para_gens):
            if self.do_mutation: self.mutate()
            self.genDrift()
    
    def genDrift(self):
        '''
        Simulates genetic drift & selection.
        '''
        if self.is_dip and not self.is_mixed_vec: return
        if not self.is_mixed: return
        all_freq = self.getAlleleFreqs()
        a = self.main_all_char
        b = self.scnd_all_char
        all_prop = all_freq/self.num_genes
        asb = self.all_sel_bias[a]
        w_avg = asb*all_prop + (1. - all_prop)
        all_prop *= asb/w_avg
        if self.is_dip: all_prop = self.rng.binomial(self.num_genes, all_prop)/self.num_genes
        if not all_prop or all_prop == 1:
            all_prop_bool = bool(all_prop)
            self.genotype_freqs[self.ploidy*a] = all_prop_bool*self.pc
            self.genotype_freqs[self.ploidy*b] = (not all_prop_bool)*self.pc
            if self.is_dip: self.genotype_freqs[a + b] = 0
            return
        if self.is_hap:
            gtfs_big = round(all_prop*self.pc_flt)
            gtfs_sml = round((1.-all_prop)*self.pc_flt)
            self.genotype_freqs[a] = gtfs_big
            self.genotype_freqs[b] = gtfs_sml
        else:
            inv_prop = 1. - all_prop
            probs = [inv_prop**2, 2*all_prop*inv_prop, all_prop**2]
            all_dist = self.rng.multinomial(n=self.pc, pvals=probs)
            self.genotype_freqs[b + b] = all_dist[0]
            self.genotype_freqs[a + b] = all_dist[1]
            self.genotype_freqs[a + a] = all_dist[2]

    def mutate(self):
        '''
        Effects the mutations observed over a single generation.
        '''
        if self.is_dip: return
        pre_gtfs = self.genotype_freqs.copy()
        num_muts = 0
        for gt in self.genotype_freqs:
            if not pre_gtfs[gt]: continue
            mut_param = self.mut_chance*pre_gtfs[gt]
            if mut_param < 0.05: continue
            if mut_param > 25: num_muts = int(mut_param)
            else: num_muts = self.rng.poisson(mut_param)
            if not num_muts: continue
            self.genotype_freqs[gt] -= num_muts
            self.genotype_freqs[gt.swapcase()] += num_muts

    def getAlleleFreqs(self):
        '''
        Gets the frequencies of each allele in the genotypes present in the individual's parasites. Used primarily in `genDrift`.
        '''
        if self.is_hap: return self.genotype_freqs[self.main_all_char]
        else:
            mac = self.main_all_char
            return 2*self.genotype_freqs[mac+mac] + self.genotype_freqs[mac+self.scnd_all_char]

    def getGenotypes(self):
        '''
        Returns all present genotypes as a list.
        '''
        return [gt for gt in self.genotype_freqs if self.genotype_freqs[gt]]
    
    def getGenotypeTransWeights_unwgt(self):
        '''
        Returns the genotypes' transmission weights.
        '''
        return [self.genotype_freqs[gt]/self.pc_flt for gt in self.genotype_freqs]

    def doesContactTransmit(self):
        '''
        Whether or not a contact between this individual and a susceptible individual resutls in a transmission.
        '''
        gtf_wgt_sum = 0
        for gt in self.genotype_freqs: gtf_wgt_sum += self.all_transm_probs[gt[0]]*self.genotype_freqs[gt]/self.pc_flt
        return random() < gtf_wgt_sum

    def infectMix(self, pc_num: int=1, do_test_contact: bool=True):
        '''
        Performs a mixed infection (dict of parasite genotype to count). Can also be used to simply get a distribution of parasites if desired.

        ### Parameters
        - `pc_num`: The number of parasites to use in the distribution.
        - `do_test_contact`: Whether or not to allow the transmission to fail, depending on the alleles' transmission probabilities.
        '''
        if do_test_contact:
            if not self.doesContactTransmit(): return dict.fromkeys(self.genotype_freqs, 0)
        return dictify(self.genotype_freqs.keys(), self.rng.multinomial(pc_num, self.getGenotypeTransWeights_unwgt()))
    
    def infectSelf(self, pc_num: int, strn: str) -> list[str]:
        '''
        Infects the individual itself.

        ### Parameters
        - `pc_num`: The number of parasites in the transmission.
        - `strn`: The strain to infect it with.
        '''
        if pc_num > self.pc: pc_num = self.pc
        to_replace = self.infectMix(pc_num, do_test_contact=False)
        strn_match = self.match(strn)
        for stn in to_replace:
            self.genotype_freqs[stn] -= to_replace[stn]
            if self.genotype_freqs[stn] < 0: to_replace[stn] += self.genotype_freqs[stn]; self.genotype_freqs[stn] = 0
            self.genotype_freqs[strn_match] += to_replace[stn]
    
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
        if self.is_dip: self.is_mixed_vec = True
        self.genotype_freqs = dict.fromkeys(self.genotype_freqs, 0)
        for strn in mix:
            matched = self.match(strn)
            amt_raw = self.pc_flt*(mix[strn]/pc_transmitted) + rem
            amt = int(amt_raw)
            rem = amt_raw - amt
            self.genotype_freqs[matched] += amt
        if rem > 0.999: self.genotype_freqs[matched] += 1

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
        else: return 1/self.num_gnts_present
    
    def correction(self, sn: str=''):
        '''
        Returns either true or false with a probability of 1/(num. strains present) for no parameter, or the frequency of the given strain.
        '''
        return random() < self.correction_det(sn)
    
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
        if self.is_hap:
            if not self.genotype_freqs[self.main_all_char]: return False
            if self.genotype_freqs[self.main_all_char] == self.pc: return False
            return True
        if not self.is_mixed_vec: return False
        been_a_something_already = False
        for gt in self.genotype_freqs:
            if been_a_something_already and self.genotype_freqs[gt]: return True
            elif self.genotype_freqs[gt]: been_a_something_already = True; continue
        return False
    
    @property
    def num_gnts_present(self):
        if self.is_hap:
            if not self.genotype_freqs[self.main_all_char]: return 1
            if self.genotype_freqs[self.main_all_char] == self.pc: return 1
            return 2
        else:
            been_a_nothing_already = False
            for gt in self.genotype_freqs:
                if been_a_nothing_already and not self.genotype_freqs[gt]: return 1
                elif not self.genotype_freqs[gt]: been_a_nothing_already = True; continue
            if been_a_nothing_already: return 2
            return 3
    
    def __str__(self):
        return f'{self.genotype_freqs}'
    
    def __repr__(self):
        return self.__str__()