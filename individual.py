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
    
    def genDrift(self, a: str):
        # a = self.allele_2_str[a_a][0]
        # ps = self.trans_ps[self.allele_freqs[a]]
        # new_num = random.choices(list(range(self.pc+1)), ps)[0]
        # self.allele_freqs[a] = new_num
        # self.allele_freqs[invChar(a)] = self.pc - new_num

        return
    
    def mutate(self):
        num_muts = self.rng.binomial(self.pc, self.mut_chance)
        for i in range(num_muts):
            mut_src = self.infect()
            pot_tgts = list(self.genotype_freqs.keys())
            pot_tgts.pop(pot_tgts.index(mut_src))
            mut_tgt = random.choice(pot_tgts)
            self.genotype_freqs[mut_src] -= 1
            self.genotype_freqs[mut_tgt] += 1
            if self.file is not None: self.file.write('\tmut\n')
    
    def reproduce(self, a: allele): # (currently) deprecated
        genes = self.allele_2_str[a]
        new_gen: dict[str, int] = {}
        for g in genes: new_gen[g] = 0
        for i in range(self.pc):
            parents = [self.getWeightedAllele(a) for i in range(2)]
            if self.is_hap: new_gen[random.choice(parents)] += 1
            else:
                gene_list = [random.choice(p) for p in parents]
                gene_list.sort()
                new_gen[''.join(gene_list)] += 1
        for g in new_gen: self.allele_freqs[g] = new_gen[g]

    def getWeightedAllele(self, a: allele): # deprecated
        return random.choices(self.allele_2_str[a], self.getAlleleDist(a))[0]
    
    def getAlleleDist(self, a: allele): # deprecated
        genes = self.allele_2_str[a]
        return [self.allele_freqs[g]/self.pc for g in genes]
    
    def getGenotypes(self):
        return [gt for gt in self.genotype_freqs if self.genotype_freqs[gt]]

    def infectMult(self, num: int=1) -> dict[str, int]:
        gtf_vals = np.array(list(self.genotype_freqs.values()))
        return dictify(self.genotype_freqs.keys(), self.rng.multinomial(num, normalise(gtf_vals)))

    def infect(self):
        return random.choices(list(self.genotype_freqs.keys()), [gtf/self.pc for gtf in list(self.genotype_freqs.values())])[0]
    
    def getAlleleFreqs(self):
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
        if self.file is not None: self.file.write('\tinfectSelf\n')
        # print('infectSelf')
        if pc_num > self.pc: pc_num = self.pc
        to_replace = self.infectMult(pc_num)
        for stn in to_replace:
            self.genotype_freqs[stn] -= to_replace[stn]
            self.genotype_freqs[self.match(strn)] += to_replace[stn]
    
    def infectSelfMult(self, mix: dict[str, int]):
        if self.file is not None: self.file.write('\tinfectSelfMult\n')
        # print('infectSelfMult')
        if sum(mix.values()) >= self.pc:
            self.setToMix(mix)
            return
        for strn in mix: self.infectSelf(mix[strn], strn)
    
    def setToMix(self, mix: dict[str, int]):
        if self.file is not None: self.file.write('\tsetToMix\n')
        # print('setToMix')
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

    def match(self, s2m: str): # matches the given strain to the format required of this individual type
        m2u = dipify
        if self.is_hap: m2u = hapify
        return m2u(s2m)
    
    def correction(self, sn: str=''):
        if sn: return random.random() < self.genotype_freqs[sn]/self.pc
        return random.random() < 1/len(self.getGenotypes())

    def storeData(self):
        if self.file is None and random.random() <= self.store_chance:
            self.file = open(f'{int(random.random()*1e6)}.dat', 'x')
            [self.file.write(f'{gnt}\t') for gnt in self.genotype_freqs]
            self.file.write('\n')
    
    def rebuild(self): # (currently) deprecated
        new_indv = individual(self.pc, tps=self.trans_ps)
        new_indv.file = self.file
        for a in self.allele_freqs: new_indv.allele_freqs[a] = self.allele_freqs[a]
        return new_indv
    
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