from gen_funcs import *
from func_lib import *
from allele import *
import random

class individual:
    '''
    The class for explicitly-modelled individuals.
    '''
    pc = 0
    allele_freqs: dict[str, int] = {}
    allele_2_str: dict[allele, list[str]] = {}
    alleles: list[allele] = []
    trans_ps = [[]]
    file = None
    is_hap = False
    do_sr = False
    mut_chance = 0.0
    para_gens = 1

    def __init__(self, alleles: list[allele]=[], gnt: str='', gdm=wf, tps: list[list[float]]=None, **kwargs):
        self.__dict__.update(kwargs)
        self.allele_freqs: dict[str, int] = {}
        self.allele_2_str: dict[allele, list[str]] = {}
        self.alleles = alleles
        for a in alleles:
            self.allele_2_str[a] = genGenotypes([a], self.is_hap)
            for a_s in self.allele_2_str[a]: self.allele_freqs[a_s] = 0
        for g in gnt.split('.'): self.allele_freqs[g] = self.pc
        if tps is None: self.trans_ps = gdm(self.pc, self.is_hap)
        else: self.trans_ps = tps

    def simPara(self):
        for i in range(self.para_gens):
            if self.file is not None:
                self.file.write('\t'.join(['\t'.join([str(int(round(self.pc*s))) for s in self.getAlleleDist(a)])
                                             for a in self.allele_2_str])+'\n')
            for a in self.alleles:
                if self.is_hap: self.genDrift(a)
                if self.mut_chance: self.mutate(a)
                if self.do_sr and not i: self.reproduce(a)
    
    def genDrift(self, a_a: allele):
        if self.is_dip: return
        a = self.allele_2_str[a_a][0]
        ps = self.trans_ps[self.allele_freqs[a]]
        new_num = random.choices(list(range(self.pc+1)), ps)[0]
        self.allele_freqs[a] = new_num
        self.allele_freqs[invChar(a)] = self.pc - new_num
    
    def mutate(self, a: allele):
        if random.random() <= self.mut_chance:
            genes = self.allele_2_str[a].copy()
            mut_src = random.choices(genes, self.getAlleleDist(a))[0]
            genes.pop(genes.index(mut_src))
            mut_tgt = random.choice(genes)
            self.allele_freqs[mut_src] -= 1
            self.allele_freqs[mut_tgt] += 1
            if self.file is None and random.random() <= self.mut_chance/20:
                self.file = open(f'{int(random.random()*1e6)}.dat', 'x')
                [self.file.write(f'{a}\t') for a in self.allele_freqs]
                self.file.write('\n')
            if self.file is not None: self.file.write('\tmut\t\n')
    
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
        return '.'.join([self.getWeightedAllele(a) for a in self.allele_2_str])
    
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