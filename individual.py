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
        gnts = genGenotypes(alleles, self.is_hap)
        for g in gnts: self.genotype_freqs[g] = 0
        if gnt: self.genotype_freqs[gnt] = self.pc
        if tps is None: self.trans_ps = gdm(self.num_genes)
        else: self.trans_ps = tps
        self.rng = np.random.default_rng()
        self.gene_range: list[int] = list(range(self.num_genes+1))
        if self.pc_to_transmit > self.pc: self.pc_to_transmit = self.pc
        self.marked_for_death = False

    def simPara(self, times: list):
        '''
        Simulates the parasites' behavior. (`times` is for speed-recording purposes; refer to sim_lib for details.)
        '''
        tba = ''
        for i in range(self.para_gens):
            if self.file is not None: self.file.write('\t'.join([str(self.genotype_freqs[g]) for g in self.genotype_freqs])+'\n')
            tm = time.time()
            allele_freqs = self.getAlleleFreqs()
            times[9] += time.time() - tm # kind of expensive
            new_allele_freqs = allele_freqs.copy()
            new_genotypes = ['']*self.pc
            curr_ind = len(allele_freqs)
            for a in allele_freqs:
                curr_ind -= 1
                tm = time.time()
                new_allele_freqs[a] = random.choices(self.gene_range, self.trans_ps[allele_freqs[a]])[0]
                times[6] += time.time() - tm
                tm = time.time()
                all_prop = new_allele_freqs[a]/self.num_genes
                probs = [1-all_prop, all_prop]
                if self.is_dip: probs = [(1-all_prop)**2, 2*all_prop*(1-all_prop), all_prop**2]
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
            if self.do_mutation: self.mutate()
            times[13] += time.time() - tm
            # reproduction placeholder
            self.storeData()
        return times
    
    def mutate(self):
        '''
        Effects the mutations observed over a single generation.
        '''
        num_muts = self.rng.binomial(self.pc, self.mut_chance)
        for i in range(num_muts):
            mut_src = self.infect()
            pot_tgts = list(self.genotype_freqs.keys())
            pot_tgts.pop(pot_tgts.index(mut_src))
            mut_tgt = random.choice(pot_tgts)
            self.genotype_freqs[mut_src] -= 1
            self.genotype_freqs[mut_tgt] += 1
            if self.file is not None: self.file.write('\tmut\n')
    
    def reproduce(self, a: allele):
        '''
        Currently deprecated; may eventually be used for sexual reproduction modelling.
        '''
        return
    
    def getGenotypes(self):
        '''
        Returns all present genotypes as a list.
        '''
        return [gt for gt in self.genotype_freqs if self.genotype_freqs[gt]]

    def infectMult(self, num: int=1) -> dict[str, int]:
        '''
        Performs multiple infections. Returns a dict of strain to # times infected.
        '''
        gtf_vals = np.array(list(self.genotype_freqs.values()))
        return dictify(self.genotype_freqs.keys(), self.rng.multinomial(num, normalise(gtf_vals)))

    def infect(self):
        '''
        Performs an infection.
        '''
        return random.choices(list(self.genotype_freqs.keys()), [gtf/self.pc for gtf in list(self.genotype_freqs.values())])[0]
    
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
        to_replace = self.infectMult(pc_num)
        for stn in to_replace:
            self.genotype_freqs[stn] -= to_replace[stn]
            self.genotype_freqs[self.match(strn)] += to_replace[stn]
    
    def infectSelfMult(self, mix: dict[str, int]):
        '''
        Infects the individual with the given strain distribution.
        '''
        if self.file is not None: self.file.write('\tinfectSelfMult\n')
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
        for strn in mix:
            matched = self.match(strn)
            amt_raw = self.pc*(mix[strn]/pc_transmitted) + rem
            amt = int(amt_raw)
            rem = amt_raw - amt
            self.genotype_freqs[matched] = amt
        if rem > 0.999: self.genotype_freqs[matched] += 1

    def match(self, s2m: str):
        '''
        Matches the given strain to the format required of this individual type (i.e., haploid or diploid).
        '''
        m2u = dipify
        if self.is_hap: m2u = hapify
        return m2u(s2m)
    
    def correction(self, sn: str=''):
        '''
        'Corrects' for overcounting by rolling a die, with success weighted either by the given strain's frequency or (if blank)
        the number of strains present.
        '''
        if sn: return random.random() < self.genotype_freqs[sn]/self.pc
        return random.random() < 1/len(self.getGenotypes())

    def storeData(self):
        '''
        Stores genotype frequency data, maybe (depends on `store_chance`).
        '''
        if self.file is None and random.random() <= self.store_chance:
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