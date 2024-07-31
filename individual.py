from gen_drift import *
import random

class individual:
    '''
    The class for explicitly-modelled individuals. Note: currently haploid only!
    '''
    pc = 0
    allele_freqs = {}
    trans_ps = [[]]
    file = None

    def __init__(self, pc: int, alleles: list[str], gdm=wf):
        self.pc = pc
        for a in alleles: self.allele_freqs[a.lower()] = pc*(a == a.lower()) # records recessive explicitly, dominant implicitly
        self.trans_ps = gdm(pc)
    
    def genDrift(self, num_gens: int=1, mut_chance: float=0.0):
        for i in range(num_gens):
            if self.file is not None:
                self.file.write('\t\t'.join(['\t\t'.join([str(int(self.pc*s)) for s in self.getAlleleDist(a)])
                                             for a in self.allele_freqs])+'\n')
            for a in self.allele_freqs:
                a_num = self.allele_freqs[a]
                ps = self.trans_ps[a_num]
                new_num = random.choices(list(range(self.pc+1)), ps)[0]
                if random.random() <= mut_chance:
                    mut_q = random.choices([-1, 1], self.getAlleleDist(a))[0]
                    new_num += mut_q
                    if self.file is None and random.random() <= 0.05:
                        self.file = open(f'{hex(id(self))}.dat', 'x')
                        [self.file.write(f'{a.lower()}\t\t{a.upper()}\n') for a in self.allele_freqs]
                self.allele_freqs[a] = new_num
    
    def getAlleleDist(self, a: str):
        return [self.allele_freqs[a]/self.pc, 1-self.allele_freqs[a]/self.pc]

    def infect(self):
        return '.'.join([random.choices([a.lower(), a.upper()], self.getAlleleDist(a))[0] for a in self.allele_freqs])