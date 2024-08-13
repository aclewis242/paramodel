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
    mut_fac = 0
    rng: np.random.Generator = None

    def __init__(self, alleles: list[allele]=[], gnt: str='', gdm=wf, tps: list[list[float]]=None, **kwargs):
        self.__dict__.update(kwargs)
        self.genotype_freqs: dict[str, int] = {}
        gnts = genGenotypes(alleles, self.is_hap)
        for g in gnts: self.genotype_freqs[g] = 0
        self.genotype_freqs[gnt] = self.pc
        if not self.do_mutation: self.mut_chance *= self.mut_fac
        if tps is None: self.trans_ps = gdm(self.num_genes)
        else: self.trans_ps = tps
        self.rng = np.random.default_rng()

    def simPara(self, times: list):
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
                new_allele_freqs[a] = random.choices(list(range(self.num_genes+1)), self.trans_ps[allele_freqs[a]])[0]
                all_prop = new_allele_freqs[a]/self.num_genes
                probs = [1-all_prop, all_prop]
                if self.is_dip: probs = [(1-all_prop)**2, all_prop*(1-all_prop), all_prop**2]
                tm = time.time()
                all_dist = self.rng.multinomial(n=self.pc, pvals=probs)
                times[10] += time.time() - tm # equally kind of expensive
                j = 0
                tm = time.time()
                for a_d_i in range(self.ploidy+1): # essentially just a pc-length for loop
                    # potential speed improvement: pre-define range (ploidy doesn't change)
                    a_d = all_dist[a_d_i]
                    for i in range(a_d):
                        idx = i + j
                        if a_d_i == 0: new_genotypes[idx] += self.ploidy*a.lower() # as for loop (/if statement?)
                            # new_genotypes[idx] += a.lower()
                            # if self.is_dip: new_genotypes[idx] += a.lower()
                        elif a_d_i == 1:
                            if self.is_hap: new_genotypes[idx] += a.upper()
                            else: new_genotypes[idx] += (a.upper() + a.lower()) # f-string?
                        elif a_d_i == 2: new_genotypes[idx] += 2*a.upper() # for/if
                        if curr_ind: new_genotypes[idx] += '.'
                    j += a_d
                times[11] += time.time() - tm
                tm = time.time()
                random.shuffle(new_genotypes)
                times[12] += time.time() - tm # most expensive but not by a huge amount
            self.genotype_freqs = self.genotype_freqs.fromkeys(self.genotype_freqs, 0)
            tm = time.time()
            for n_g in new_genotypes: self.genotype_freqs[n_g] += 1
            times[12] += time.time() - tm
            tm = time.time()
            self.mutate()
            times[13] += time.time() - tm
            # reproduction placeholder
        return times
    
    def genDrift(self, a: str):
        # a = self.allele_2_str[a_a][0]
        # ps = self.trans_ps[self.allele_freqs[a]]
        # new_num = random.choices(list(range(self.pc+1)), ps)[0]
        # self.allele_freqs[a] = new_num
        # self.allele_freqs[invChar(a)] = self.pc - new_num

        return
    
    def mutate(self):
        if random.random() <= self.mut_chance:
            mut_src = self.infect()
            pot_tgts = list(self.genotype_freqs.keys())
            pot_tgts.pop(pot_tgts.index(mut_src))
            mut_tgt = random.choice(pot_tgts)
            self.genotype_freqs[mut_src] -= 1
            self.genotype_freqs[mut_tgt] += 1
            self.storeData()
            if self.file is not None: self.file.write('\tmut\n')
    
    def reproduce(self, a: allele):
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

    def getWeightedAllele(self, a: allele):
        return random.choices(self.allele_2_str[a], self.getAlleleDist(a))[0]
    
    def getAlleleDist(self, a: allele):
        genes = self.allele_2_str[a]
        return [self.allele_freqs[g]/self.pc for g in genes]

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
    
    def storeData(self):
        mf = self.mut_fac
        if self.do_mutation: mf /= self.mut_fac
        if self.file is None and random.random() <= self.mut_chance/(200*mf):
            self.file = open(f'{int(random.random()*1e6)}.dat', 'x')
            [self.file.write(f'{gnt}\t') for gnt in self.genotype_freqs]
            self.file.write('\n')
    
    def rebuild(self):
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