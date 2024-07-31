from func_lib import *
from individual import *

class population:
    '''
    The class for a population.
    '''
    sus = -1
    inf = {}
    rec = {}
    pn = ''
    tot_pop = 0
    is_vector = False
    indvs = []
    pc = 0

    def __init__(self, p0: list[int], pn: str='', isn: str='init', pc: int=120):
        '''
        Initialises the population.

        ### Parameters
        - `p0`: A 3-element integer list of the starting population amounts (S, I, R)
        - `pn`: The population's name
        - `isn`: The name of the initial strain
        '''
        self.sus = p0[0] + 1 # adding 1 to ensure the population is always non-zero (irrelevant in the grand scheme of things)
        self.inf = {}
        self.rec = {}
        self.inf[isn] = p0[1]
        self.rec[isn] = p0[2]
        self.pn = pn
        self.tot_pop = sum(p0)
        self.pc = pc
    
    def getPop(self, sn: str='init') -> list[int]:
        '''
        Returns a 3-element S, I, R list for the given strain.
        '''
        return [self.sus, self.inf[sn], self.rec[sn]]
    
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
        self.sus += p[0]
        self.inf[sn] += p[1]
        if self.is_vector:
            if p[1] >= 0:
                for i in range(int(p[1])): self.indvs += [individual(self.pc, sn.split('.'))]
            else: self.individuals = self.indvs[:int(p[1])]
            random.shuffle(self.indvs)
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