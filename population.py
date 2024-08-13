from func_lib import *
from individual import *

class population:
    '''
    The class for a population.
    '''
    sus = -1
    inf: dict[str, int] = {}
    rec: dict[str, int] = {}
    pn = ''
    tot_pop = 0
    is_vector = False
    indvs: list[individual] = []
    is_hap = False
    do_indvs = False
    indv_params = {} # pc, is_hap, do_sr, mut_chance, alleles

    def __init__(self, p0: list[int], pn: str='', isn: str='init', **kwargs):
        '''
        Initialises the population.

        ### Parameters
        - `p0`: A 3-element integer list of the starting population amounts (S, I, R)
        - `pn`: The population's name
        - `isn`: The name of the initial strain
        '''
        self.sus = p0[0] + 1 # adding 1 to ensure the population is always non-zero (irrelevant in the grand scheme of things)
        self.inf: dict[str, int] = {}
        self.rec: dict[str, int] = {}
        self.indvs: list[individual] = []
        self.inf[isn] = p0[1]
        self.rec[isn] = p0[2]
        self.pn = pn
        self.tot_pop = sum(p0)
        self.__dict__.update(kwargs)
        self.indv_params = kwargs
    
    def getPop(self, sn: str='init') -> list[int]:
        '''
        Returns a 3-element S, I, R list for the given strain.
        '''
        I = self.inf[sn]
        return [self.sus+sum(self.inf.values()), self.inf[sn], self.rec[sn]]
    
    def getAllPop(self):
        '''
        Returns all population elements as a list. S is first, then all Is, then all Rs.
        '''
        return [self.sus] + list(self.inf.values()) + list(self.rec.values())
    
    def getAllPopNms(self):
        '''
        Returns the names of all population elements. S is first, then all Is, then all Rs.
        '''
        return ['S'] + [f'I ({sn})' for sn in self.inf.keys()] + [f'R ({sn})' for sn in self.rec.keys()]
    
    def addPop(self, p: list[int], sn: str='init'):
        '''
        Adds the given population quantities.

        ### Parameters
        - `p`: The quantities to add, as a 3-element list (S, I, R).
        - `sn`: The strain to add them to.
        '''
        sn = self.match(sn)
        self.sus += p[0]
        if self.inf[sn] + p[1] < 0: p[1] = -self.inf[sn]
        self.inf[sn] += p[1]
        if self.do_indvs:
            if p[1] >= 0:
                for i in range(int(p[1])): self.indvs += [individual(gnt=sn, **self.indv_params)]
            else: self.indvs = shuffle(self.indvs[:int(p[1])])
        if self.rec[sn] + p[2] < 0: p[2] = -self.rec[sn]
        self.rec[sn] += p[2]
        self.tot_pop += (p[0] + p[1] + p[2])
    
    @property
    def individuals(self) -> list[individual]:
        return self.indvs
    
    @individuals.setter
    def individuals(self, value: list[individual]):
        self.indvs = value

    def addStrain(self, nsn: str):
        '''
        Creates a new strain with the given name.
        '''
        if nsn in self.inf.keys(): return
        self.inf[nsn] = 0
        self.rec[nsn] = 0

    @property
    def is_dip(self):
        return not self.is_hap
    
    @is_dip.setter
    def is_dip(self, value: bool):
        self.is_hap = not value
    
    def match(self, s2m: str): # matches the given strain to the format required of this population
        m2u = dipify
        if self.is_hap: m2u = hapify
        return m2u(s2m)

    def printDat(self):
        '''
        Prints the object's information to the console.
        '''
        print(''.join([f'Population {self.pn}\n',
                        f'S:\t{float2SN(self.sus)}\n',
                        f'\tI\n',
                        '\n'.join([f'{k}:\t{float2SN(self.inf[k])}' for k in self.inf]),
                        f'\n\tR\n',
                        '\n'.join([f'{k}:\t{float2SN(self.rec[k])}' for k in self.rec])]))

    def __str__(self):
        return self.pn
    
    def __repr__(self):
        return self.__str__()