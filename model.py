from func_lib import *
from population import *
from allele import *
import random

class SIR:
    '''
    The stochastic SIR model class.
    '''
    pn = ''
    pn_full = ''
    sn = 'init'
    genotype = ''
    alleles = []
    bd = -1.0
    ir = -1.0
    rr = -1.0
    wi = -1.0
    itr: dict[population, float] = {}
    itr_keys: list[population] = []
    mr = 0.0
    dt = 0.0
    bds = 0.0
    Rs = []
    Es = []
    num_Es = -1
    is_vector = False
    do_mixed_infs = False
    
    def __init__(self, p0: population, **kwargs):
        '''
        Initialises the model with the given parameters.

        ### Parameters
        - `p0`: Initial population, as a 3-element list (S, I, R)
        - `pn`: The name (short) of the population this model belongs to
        - `pn_full`: The full name of the population this model belongs to
        - `sn`: The name of the strain this model belongs to (i.e. genotype, in most cases)
        - `bd`: Birth/death rate
        - `ir`: Infection rate
        - `rr`: Recovery rate
        - `wi`: Waning immunity rate
        - `itr`: Interspecific transmission rates from this model to other populations (dict)
        - `mr`: Mutation rate
        - `is_vector`: Whether or not this model describes a vector population (bool)
        - `do_mixed_infs`: Whether or not this model does mixed infections
        '''
        self.pop = p0
        self.__dict__.update(kwargs)
        self.pop.is_vector = self.is_vector
        self.pop.pn = self.pn
        self.pop.do_mixed_infs = self.do_mixed_infs
        self.pop.addStrain(self.sn)
        E1 = [1, 0, 0]  # Birth (used 3 times, one for each pop type)
        E2 = [-1, 1, 0] # Infection
        E3 = [0, -1, 1] # Recovery
        E4 = [-1, 0, 0] # Death of susceptible
        E5 = [0, -1, 0] # Death of infected
        E6 = [0, 0, -1] # Death of recovered
        E7 = [1, 0, -1] # Waning immunity
        self.Es = [E1, E2, E3, E4, E5, E6, E7, E1, E1, E1] # first E1 is 'legacy'
        self.num_Es = len(self.Es)
        self.itr_keys = list(self.itr.keys())
        self.setRs()

    def setRs(self):
        '''
        Generates the different transition rates based on the model's parameters and population.
        '''
        # [S, I, R] = self.pop.getPop(self.sn)
        # S = float(int(S))
        # I = float(int(I))
        # R = float(int(R))
        # N = S + I + R
        # if not N: return [0. for r in self.Rs]
        S = self.pop.sus
        R = self.pop.rec[self.sn]
        I_UW = 0
        I_WS = 0
        inf_indvs = self.pop.getSusInf(self.sn, is_present=True)
        for ind in inf_indvs:
            I_UW += ind.correction_det()            # unweighted infections (simple 'yes/no' on strain presence)
            I_WS += ind.correction_det(sn=self.sn)  # weighted according to how prevalent the strain is inside the indv
        self.Rs = [0, # N*self.bd*0,     # 'legacy', since births are now pop-delineated (like deaths always were)
                    0, # self.ir*S*I/N,  # also 'legacy' (infections are now handled by vectors, not intra-pop things)
                    self.rr*I_WS,
                    self.bds*S,
                    self.bd*I_UW,
                    self.bd*R,
                    self.wi*R,
                    self.bds*S,
                    self.bd*I_UW,
                    self.bd*R] + [self.itr[p2]*I_WS*(p2.sus+p2.getSusInfNum(self.sn))/(self.pop.tot_pop+p2.tot_pop) for p2 in self.itr]
                    # is the tot_pop # too high?
        return self.Rs

    def trans(self, idx: int, rpt: int=1):
        '''
        Effects the changes in the population dictated by the simulation.

        ### Parameters
        - `idx`: The index of the desired event, corresponding with the order of `Rs`.
        - `rpt`: The number of times to repeat said event.
        '''
        pop = self.pop
        pc_trans_src = self.pop.pc_to_transmit
        if idx >= self.num_Es:
            pop = self.itr_keys[idx-self.num_Es]
            idx = 1
            num_inf = 0
            num_mixes = 0
            num_loops = 0
            max_loops = 10000
            while num_inf < rpt:
                random.shuffle(self.pop.individuals)
                for indv in self.pop.individuals:
                    if indv.correction(sn=self.sn):
                        num_inf += 1
                        if indv.is_mixed:
                            pop.infectMix(indv.infectMult(pc_trans_src))
                            num_mixes += 1
                    if num_inf == rpt: break
                num_loops += 1
                if num_loops >= max_loops: break
            rpt -= num_mixes
        if rpt: pop.addPop(list(np.multiply(self.Es[idx], rpt)), self.sn, pc_trans_src)
    
    def newStrain(self, nsn='new'):
        '''
        Generates a copy of this model with the given strain name.
        '''
        new_mdl = SIR(self.pop, sn=nsn, **self.__dict__)
        new_mdl.itr = dict(new_mdl.itr)
        new_mdl.itr_keys = list(new_mdl.itr.keys())
        return new_mdl
    
    def mutate(self, param: str, fac: float, vec: 'SIR'=None):
        '''
        Effects the given parameter change.

        ### Parameters
        - `param`: The parameter of the model to change.
        - `fac`: The numerical factor to change it by. This value is used directly.
        - `vec`: The corresponding (same strain) model for the vector population. Only necessary if `param` is `itr`.
        '''
        if type(self.__dict__[param]) is dict:
            for k in self.__dict__[param]: self.__dict__[param][k] *= fac
            if vec is not None: vec.__dict__[param][self.pop] *= fac
        else: self.__dict__[param] *= fac
    
    def mutateMult(self, params: list[str], fac: float, vec: 'SIR'=None):
        '''
        Effects the given parameter changes.

        ### Parameters
        - `params`: The parameters of the model to change.
        - `fac`: The numerical factor to change them by (all the same). This value is used directly.
        - `vec`: The corresponding (same strain) model for the vector population. Only necessary if `params` includes `itr`.
        '''
        for p in params: self.mutate(p, fac, vec)
    
    def updateGenotype(self, g: str, alleles: list[allele], vec: 'SIR'=None):
        '''
        Generates a new model based on the given genotype.

        ### Parameters
        - `g`: The genotype, as a string of characters corresponding to alleles.
        - `alleles`: The list of all possible alleles, as allele objects.
        - `vec`: The corresponding (same strain) model for the vector population. Only necessary if one of the alleles affects `itr`.
        '''
        new_model = self.newStrain(g)
        new_model.genotype = g
        for a in alleles:
            if a.char in g:
                if a.fav_pop == new_model.pn: new_model.mutate(a.param, 1+a.fac, vec)
                if a.unf_pop == new_model.pn: new_model.mutate(a.param, 1/(1+a.fac), vec)
        return new_model

    def r0(self, vec_mdl: 'SIR') -> float:
        '''
        Estimates R0 for the given model. Meant to be more a vague guideline than a hard and fast rule. Takes the corresponding (same strain)
        vector model as an input.
        '''
        return (vec_mdl.pop.tot_pop/self.pop.tot_pop)*self.itr[vec_mdl.pop]*vec_mdl.itr[self.pop]/(self.rr*vec_mdl.bd)

    def printParams(self):
        '''
        Prints each of the model's parameters and their values to the console.
        '''
        print(f'\n{self.__str__()} params:')
        [print(f'{k}: {self.__dict__[k]}') for k in self.__dict__]
    
    def __str__(self):
        return f'population {self.pn}, strain {self.sn}'
    
    def __repr__(self):
        return self.__str__()