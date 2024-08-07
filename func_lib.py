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
    if f < 1e3: return f'{f}'
    else:
        pwr = np.floor(np.log10(f))
        return f'{int(f/(10**(pwr-p)))/(10**p)}e{int(pwr)}'

def printMat(m: list[list]):
    rv = '\n'.join(['\t'.join([str(int(100*i)) for i in j]) for j in m])
    print(rv)
    return rv

def shuffle(l: list):
    order = list(range(len(l)))
    random.shuffle(order)
    return [l[i] for i in order]