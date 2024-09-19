import random
import numpy as np

def normalise(l: np.ndarray[float]):
    '''
    Normalises the given list such that its sum is 1.
    '''
    return l/sum(l)

def fn(f_n: str):
    '''
    Converts the given string into a proper file name.
    '''
    return f_n.lower().replace(' ', '_')

def float2SN(f: float, p: int=2):
    '''
    Converts the given number to scientific notation with the given precision.

    ### Parameters
    - `f`: The number in question.
    - `p`: The number (int) of digits of precision (past the decimal point). For example: f=1302, p=2 -> 1.30e3.
    '''
    if f < 10**(p+1): return f'{f}'
    else:
        pwr = np.floor(np.log10(f))
        return f'{int(f/(10**(pwr-p)))/(10**p)}e{int(pwr)}'

def roundNum(f: float, prec: int=2):
    return round(f*(10**prec))/(10**prec)

def normPercentList(l: list[float]) -> list[float]:
    return list(100*normalise(np.array(l)))

def printFloatList(l: list[float]):
    [print(f'{i}:\t{roundNum(l[i])}') for i in range(len(l))]

def printMat(m: list[list]):
    '''
    Prints the given 2d matrix (usually of transition probabilities, though not necessarily) in an easier-to-read fashion.
    '''
    rv = '\n'.join(['\t'.join([str(int(100*i)) for i in j]) for j in m])
    print(rv)
    return rv

def hapify(g: str):
    '''
    Turns the given genotype into its haploid equivalent (e.g., AA.Bb.cc becomes A.B.c). Haploid genotypes are unaffected.
    '''
    return '.'.join([s[0] for s in g.split('.')])

def dipify(g: str):
    '''
    Turns the given genotype into its diploid equivalent (e.g., A.B.c becomes AA.BB.cc). Diploid genotypes are unaffected.
    '''
    return '.'.join([''.join(2*[s])[:2] for s in g.split('.')])

def listify(a) -> list:
    '''
    Produces a deep copy of type list from a given multi-dimensional list (or NumPy array). If the elements of the list are complex
    types, then these objects will still be the same in memory (the deep copy only extends to the lists themselves).
    '''
    if (type(a) is np.ndarray) or (type(a) is list): return [listify(el) for el in a]
    else: return a

def dictify(ks: list, vs: list):
    '''
    Turns the given lists of keys and values into a dict. These should be the same length!
    '''
    return {k: v for k, v in zip(ks, vs)}