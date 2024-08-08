from gen_drift import *
from func_lib import *
import random

class individual:
    '''
    The class for explicitly-modelled individuals. Note: currently haploid only!
    '''
    pc = 0
    allele_freqs: dict[str, int] = {}
    trans_ps = [[]]
    file = None
    is_hap = False
    do_sr = False
    mut_chance = 0.0
    para_gens = 1

    def __init__(self, alleles: list[str]=[], gdm=wf, tps: list[list[float]]=None, **kwargs):
        self.__dict__.update(kwargs)
        self.allele_freqs: dict[str, int] = {}
        for a in alleles: self.allele_freqs[a.lower()] = self.pc*(a == a.lower()) # records recessive explicitly, dominant implicitly
        if tps is None: self.trans_ps = gdm(self.pc, self.is_hap)
        else: self.trans_ps = tps

    def simPara(self):
        for i in range(self.para_gens):
            if self.file is not None:
                self.file.write('\t\t'.join(['\t\t'.join([str(int(round(self.pc*s))) for s in self.getAlleleDist(a)])
                                             for a in self.allele_freqs])+'\n')
            for a in self.allele_freqs:
                self.genDrift(a)
                self.mutate(a)
    
    def genDrift(self, a: str):
        a_num = self.allele_freqs[a]
        ps = self.trans_ps[a_num]
        new_num = a_num
        if self.is_hap: new_num = random.choices(list(range(self.pc+1)), ps)[0]
        self.allele_freqs[a] = new_num
    
    def mutate(self, a: str):
        if random.random() <= self.mut_chance:
            self.allele_freqs[a] += random.choices([-1, 1], self.getAlleleDist(a))[0]
            if self.file is None and random.random() <= self.mut_chance/20:
                self.file = open(f'{int(random.random()*1e6)}.dat', 'x')
                [self.file.write(f'{a.lower()}\t\t{a.upper()}\n') for a in self.allele_freqs]
            if self.file is not None: self.file.write('\tmut\t\n')
    
    def getAlleleDist(self, a: str):
        return [self.allele_freqs[a]/self.pc, 1-self.allele_freqs[a]/self.pc]

    def infect(self):
        return '.'.join([random.choices([a.lower(), a.upper()], self.getAlleleDist(a))[0] for a in self.allele_freqs])
    
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