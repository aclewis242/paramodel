from allele import *
import math

def wf(pc: int):
    '''
    Produces the Wright-Fisher genetic drift transmission matrix for the given # of alleles.
    '''
    tp_dim = pc + 1
    trans_ps = []
    for i in range(tp_dim):
        trans_ps_row = []
        for j in range(pc+1): trans_ps_row += [0]
        trans_ps += [trans_ps_row]
    for i in range(tp_dim):
        for j in range(tp_dim):
            trans_ps[i][j] = math.comb(pc, j)*((i/pc)**j)*(1 - i/pc)**(pc-j)
    return trans_ps

def genGenotypes(alleles: list[allele], is_haploid: bool):
    '''
    Generates all possible genotypes.

    ### Parameters
    - `alleles`: The list of allele objects to generate the genotypes from.
    - `is_haploid`: Whether or not the genotypes should be haploid.
    '''
    loci = ''
    for a in alleles:
        if a.locus not in loci: loci += a.locus
    num_combs = 3 - is_haploid
    gt = ['' for i in range(num_combs**len(loci))]
    for l in loci:
        for i in range(len(gt)):
            if is_haploid:
                if i%2 == 0: gt[i] += l.upper()
                if i%2 == 1: gt[i] += l.lower()
            else:
                if i%3 == 0: gt[i] += (l.upper() + l.lower())
                if i%3 == 1: gt[i] += 2*l.upper()
                if i%3 == 2: gt[i] += 2*l.lower()
            gt[i] += '.'
        gt.sort()
    gt = [g[:-1] for g in gt]
    return gt

def genAlleles(al: str):
    is_haploid = len(al) == 1
    al = al[0]
    if is_haploid: return [al.lower(), al.upper()]
    else: return [2*al.lower(), (al.upper() + al.lower()), 2*al.upper()]