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
    mr = 0.0
    dt = 0.0
    bds = -1.0
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
        self.setRs()

    def setRs(self):
        '''
        Generates the different transition rates based on the model's parameters and population.
        '''
        [S, I, R] = self.pop.getPop(self.sn)
        S = float(int(S))
        I = float(int(I))
        R = float(int(R))
        N = S + I + R
        if not N: return [0. for r in self.Rs]
        self.bds = self.bd/len(self.pop.inf)
        self.Rs = [N*self.bd*0,     # 'legacy', since births are now pop-delineated (like deaths always were)
                    self.ir*S*I/N,
                    self.rr*I,
                    self.bds*S,
                    self.bd*I,
                    self.bd*R,
                    self.wi*R,
                    self.bds*S,
                    self.bd*I,
                    self.bd*R] + [self.itr[p2]*I*p2.sus/(N+p2.tot_pop) for p2 in self.itr]
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
        def addPopMult(idx, rpt, sn):
            pop.addPop(list(map(lambda x: float(rpt)*x, self.Es[idx])), sn, pc_trans_src)
        if idx >= self.num_Es:
            pop = list(self.itr.keys())[idx-self.num_Es]
            idx = 1
            if len(self.pop.individuals): # consider moving to a bool (do_indvs) for speed
                to_infect = {}
                for i in range(int(rpt)):
                    indv = random.choice(self.pop.individuals)
                    # potential speed increase: run infectMult earlier, feed into normal way if all same gt
                    if not (self.pop.do_mixed_infs and pop.do_mixed_infs) or len(indv.getGenotypes()) == 1:
                        if indv.correction(self.sn):
                            strn = indv.infect()
                            f = open('inf_events_raw.dat', 'a')
                            f.write(f'gtfs {indv.genotype_freqs} -> {strn}\n')
                            f.close()
                            if strn in to_infect.keys(): to_infect[strn] += 1
                            else: to_infect[strn] = 1
                    elif indv.correction(self.sn): pop.infectMix(indv.infectMult(indv.pc_to_transmit))
                for strn in to_infect: addPopMult(idx, to_infect[strn], strn)
                return self.pop.getPop(self.sn)
        elif (idx == 2 or idx == 8) and pop.do_indvs:
            # note that this solution to the unbalanced inf births/deaths problem is adequate, not perfect
            random.shuffle(pop.individuals)
            chng = 0
            for indv in pop.individuals[:int(rpt)]: chng += indv.correction(self.sn*(idx==2))
            rpt = chng
        addPopMult(idx, rpt, self.sn)
        return self.pop.getPop(self.sn)
    
    def newStrain(self, nsn='new'):
        '''
        Generates a copy of this model with the given strain name.
        '''
        new_mdl = SIR(self.pop, sn=nsn, **self.__dict__)
        new_mdl.itr = dict(new_mdl.itr)
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
            # note for future: in a more complex use-case (i.e. more 'multi-model parameters'), it may be a good idea to
            # reconsider using the non-primitive number object idea from earlier (more general)
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
        def effectMutation(al: allele, mdl: 'SIR', v_m: 'SIR'=vec):
            if al.fav_pop == mdl.pn: mdl.mutate(al.param, 1+al.fac, v_m)
            if al.unf_pop == mdl.pn: mdl.mutate(al.param, 1/(1+al.fac), v_m)
            return mdl
        for a in alleles:
            if a.char in g:
                new_model = effectMutation(a, new_model)
                if self.pop.do_sel_bias and vec is not None:
                    temp_mdl = self.newStrain(g)
                    temp_mdl.r0(vec, set_biases=True, sn=a.locus)
                    effectMutation(a, temp_mdl, v_m=None).r0(vec, set_biases=True, sn=a.char)
        return new_model

    def r0(self, vec_mdl: 'SIR', set_biases: bool=False, sn: str='') -> float:
        '''
        Estimates R0 for the given model. Meant to be more a vague guideline than a hard and fast rule. Takes the corresponding (same strain)
        vector model as an input.
        '''
        r0_val = (vec_mdl.pop.tot_pop/self.pop.tot_pop)*self.itr[vec_mdl.pop]*vec_mdl.itr[self.pop]/(self.rr*vec_mdl.bd)
        if self.pop.do_sel_bias and set_biases:
            if not sn: sn = self.sn
            self.pop.gnt_sel_bias[sn] = r0_val
            if sn in vec_mdl.pop.sel_bias_lst: vec_mdl.pop.sel_bias_lst[sn] += [r0_val]
            else: vec_mdl.pop.sel_bias_lst[sn] = [r0_val]
        return r0_val

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