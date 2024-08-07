from gen_drift import *
from func_lib import *
import random

class individual:
    '''
    The class for explicitly-modelled individuals. Note: currently haploid only!
    '''
    pc = 0
    allele_freqs = {}
    trans_ps = [[]]
    file = None

    def __init__(self, pc: int, alleles: list[str]=[], gdm=wf, tps: list[list[float]]=None):
        self.pc = pc
        self.allele_freqs = {}
        for a in alleles: self.allele_freqs[a.lower()] = pc*(a == a.lower()) # records recessive explicitly, dominant implicitly
        if tps is None: self.trans_ps = gdm(pc)
        else: self.trans_ps = tps

    def genDrift(self, num_gens: int=1, mut_chance: float=0.0):
        for i in range(num_gens):
            if self.file is not None:
                self.file.write('\t\t'.join(['\t\t'.join([str(int(round(self.pc*s))) for s in self.getAlleleDist(a)])
                                             for a in self.allele_freqs])+'\n')
            for a in self.allele_freqs:
                a_num = self.allele_freqs[a]
                ps = self.trans_ps[a_num]
                new_num = random.choices(list(range(self.pc+1)), ps)[0]
                self.allele_freqs[a] = new_num
                if abs(new_num - a_num) > self.pc-2: print(f'bruh moment: going from {a_num} {a} to {new_num}')
                if random.random() <= mut_chance:
                    self.allele_freqs[a] += random.choices([-1, 1], self.getAlleleDist(a))[0]
                    if self.file is None and random.random() <= mut_chance/20:
                        self.file = open(f'{int(random.random()*1e6)}.dat', 'x')
                        [self.file.write(f'{a.lower()}\t\t{a.upper()}\n') for a in self.allele_freqs]
                    if self.file is not None: self.file.write(f'{new_num}\t->\t{self.allele_freqs[a]}\n')
    
    def getAlleleDist(self, a: str):
        return [self.allele_freqs[a]/self.pc, 1-self.allele_freqs[a]/self.pc]

    def infect(self):
        return '.'.join([random.choices([a.lower(), a.upper()], self.getAlleleDist(a))[0] for a in self.allele_freqs])
    
    def rebuild(self):
        new_indv = individual(self.pc, tps=self.trans_ps)
        new_indv.file = self.file
        for a in self.allele_freqs: new_indv.allele_freqs[a] = self.allele_freqs[a]
        return new_indv