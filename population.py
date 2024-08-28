from func_lib import *
from individual import *
import statistics as stats

class population:
    '''
    The class for a population.
    '''
    sus = -1
    inf: dict[str, int] = {}
    rec: dict[str, int] = {}
    pn = ''
    is_vector = False
    indvs: list[individual] = []
    is_hap = False
    do_indvs = False
    indv_params = {} # pc, is_hap, do_sr, mut_chance, alleles
    rng: np.random.Generator = None
    do_mixed_infs = False
    pc_to_transmit = 0
    store_chance = 0.
    gnt_sel_bias: dict[str, float] = {}
    sel_bias_lst: dict[str, list[float]] = {}
    do_sel_bias = True

    def __init__(self, p0: list[int], pn: str='', isn: str='init', **kwargs):
        '''
        Initialises the population.

        ### Parameters
        - `p0`: A 3-element integer list of the starting population amounts (S, I, R)
        - `pn`: The population's name
        - `isn`: The name of the initial strain
        - `is_vector`: Whether or not this population is a vector
        - `is_hap`: Whether or not this population is haploid
        - `do_indvs`: Whether or not this population models infected individuals explicitly
        - `do_mixed_infs`: Whether or not this population does mixed infections
        '''
        self.sus = p0[0] + 1 # adding 1 to ensure the population is always non-zero (irrelevant in the grand scheme of things)
        self.inf: dict[str, int] = {}
        self.rec: dict[str, int] = {}
        self.indvs: list[individual] = []
        self.inf[isn] = p0[1]
        self.rec[isn] = p0[2]
        self.pn = pn
        self.__dict__.update(kwargs)
        self.indv_params = kwargs
        self.rng = np.random.default_rng()
        if not self.do_indvs: self.do_mixed_infs = False
        self.gnt_sel_bias: dict[str, float] = {}
        self.sel_bias_lst: dict[str, list[float]] = {}
    
    def getPop(self, sn: str='init') -> list[int]:
        '''
        Returns a 3-element S, I, R list for the given strain.
        '''
        return [self.sus, self.inf[sn], self.rec[sn]]
    
    def getAllPop(self, weight: bool=False):
        '''
        Returns all population elements as a list. S is first, then all Is, then all Rs.
        '''
        infs = self.inf
        if weight and self.do_indvs:
            infs = dict.fromkeys(infs, 0)
            for ind in self.indvs:
                for gt in ind.genotype_freqs: infs[gt] += ind.genotype_freqs[gt]/ind.pc
        return [self.sus] + list(infs.values()) + list(self.rec.values())
    
    def getAllPopNms(self):
        '''
        Returns the names of all population elements. S is first, then all Is, then all Rs.
        '''
        return ['S'] + [f'I ({sn})' for sn in self.inf.keys()] + [f'R ({sn})' for sn in self.rec.keys()]
    
    def addPop(self, p: list[int], sn: str='init', pc_trans: int=0):
        '''
        Adds the given population quantities.

        ### Parameters
        - `p`: The quantities to add, as a 3-element list (S, I, R).
        - `sn`: The strain to add them to.
        - `pc_trans`: The number of parasites a mixed infection uses. Note that this is based on the *source,* not the destination!
        '''
        sn = self.match(sn)
        old_data = self.printDatStr()
        test_add = [self.sus, self.inf[sn], self.rec[sn]]
        neg = -1
        for i in range(3):
            if test_add[i] < 0:
                print(test_add)
                print(f'strain {sn}')
                self.printDat()
                exit()
            if test_add[i] + p[i] < 0: neg = test_add[i]/abs(p[i])
        if neg >= 0:
            for i in range(3):
                p[i] *= neg
                if test_add[i] + p[i] < 0:
                    print(test_add)
                    print(p)
                    print(neg)
                    print(f'strain {sn}')
                    self.printDat()
                    exit()
        S = int(p[0])
        I = int(p[1])
        R = int(p[2])
        sus_new = self.sus
        inf_new = self.inf.copy()
        rec_new = self.rec.copy()
        sus_indvs: list[individual] = []
        strain_indvs: dict[str, list[individual]] = {}
        sus_pops: np.ndarray[float]
        do_shuffle = False
        if self.do_mixed_infs and (S < 0 or I):
            sus_indvs = self.getSusInf(sn)
            for sn_i in self.inf: strain_indvs[sn_i] = self.getSusInf(sn_i, is_present=True, indvs_lst=sus_indvs)
            sus_pops_lst = [self.sus]+[(sn_i != sn)*len(strain_indvs[sn_i]) for sn_i in self.inf]
            sus_pops = np.array(sus_pops_lst)/sum(sus_pops_lst)
        
        sus_new += S
        inf_new[sn] += I
        inf_indvs = []
        I_corr = 0
        if I > 0:
            if self.do_indvs:
                tba: list[individual] = []
                if not self.do_mixed_infs: tba = self.makeIndvs(sn, I)
                else: # 'mixed' in the sense of 'mixed final parasite genotypes,' not 'mixed transmission'
                    i_change, to_change = self.getChanges(I, sus_pops)
                    tba = self.makeIndvs(sn, i_change)
                    for strn in to_change:
                        indvs_to_infect = strain_indvs[strn]
                        if do_shuffle: random.shuffle(indvs_to_infect)
                        for ind in indvs_to_infect[:to_change[strn]]: ind.infectSelf(pc_trans, sn)
                        do_shuffle = True
                self.indvs += tba
        elif I < 0:
            I_corr = 0
            self.refresh()
            random.shuffle(self.indvs)
            inf_indvs = self.getSusInf(sn, is_present=True)
            for ind in inf_indvs[:-I]:
                if R or ind.correction():
                    ind.marked_for_death = True
                    other_gts = list(set(ind.getGenotypes()) - set([sn]))
                    for gt in other_gts: inf_new[gt] -= 1
                else: I_corr += 1
            self.refresh()
            inf_new[sn] += I_corr
        rec_new[sn] += R

        self.sus = sus_new
        self.inf = inf_new.copy()
        self.rec = rec_new.copy()
        self.refresh()
        f = open('last_event.dat', 'w')
        f.write(f'\n-----\nadded {p} with strain {sn}\nold data:\n')
        f.write(old_data)
        f.write('\nnew data:\n')
        f.write(self.printDatStr())
        f.close()
        # if (not self.is_vector and self.tot_pop > 1051) or len(self.indvs) < max(self.inf.values()):
        #     print(f'p: {p}, sn: {sn}')
        #     print('OLD:')
        #     print(old_data)
        #     print(f'len old indvs: {len(old_indvs)}; old sus: {old_sus}')
        #     print('NEW:')
        #     self.printDat()
        #     print(f'len new indvs: {len(self.indvs)}; new sus: {self.sus}')
        #     print(sum(self.rec.values()))
        #     [print(ind.genotype_freqs) for ind in inf_indvs]
        #     print(f'I_corr: {I_corr}')
        #     exit()
    
    def getChanges(self, pop_num: int, weights: np.ndarray[float]):
        '''
        Distributes infections between full-susceptibles and infected-susceptibles.
        '''
        to_change_lst = self.rng.multinomial(pop_num, weights)
        return to_change_lst[0], dictify(self.inf.keys(), to_change_lst[1:])
    
    def makeIndvs(self, sn: str, num_indvs: int):
        '''
        Make the given number of individuals with the given strain.
        '''
        self.store_chance /= (1+sum(self.inf.values()))
        return [individual(gnt=sn, **self.indv_params) for i in range(int(num_indvs))]

    def addStrain(self, nsn: str):
        '''
        Creates a new strain with the given name.
        '''
        if nsn in self.inf.keys(): return
        self.inf[nsn] = 0
        self.rec[nsn] = 0

    def match(self, s2m: str):
        '''
        Matches the given strain to the format required of this population (i.e., haploid or diploid).
        '''
        m2u = dipify
        if self.is_hap: m2u = hapify
        return m2u(s2m)
    
    def getSusInf(self, sn: str, is_present: bool=False, indvs_lst: list[individual]=[]):
        # get list of infected individuals either with the strain (T) or susceptible to infection by the strain (F)
        '''
        Gets list of infected individuals either with strain `sn` (`is_present`) or susceptible to infection by the strain
        (`not is_present`). A list to draw from can be provided, if desired (otherwise, it'll draw from `self.indvs`).
        '''
        is_present -= 1
        if not len(indvs_lst): indvs_lst = self.indvs
        return [ind for ind in indvs_lst if bool(ind.genotype_freqs[sn])+is_present]
    
    def getSusInfNum(self, sn: str, is_present: bool=False, indvs_lst: list[individual]=[]):
        '''
        Gets the number of infected individuals either with strain `sn` (`is_present`) or susceptible to infection by the strain
        (`not is_present`). A list to draw from can be provided, if desired (otherwise, it'll draw from `self.indvs`).
        '''
        rv = 0
        if not len(indvs_lst): indvs_lst = self.indvs
        is_present -= 1
        for ind in indvs_lst:
            if bool(ind.genotype_freqs[sn])+is_present: rv += 1
        return rv

    def infectMix(self, mix: dict[str, int]):
        '''
        Perform a mixed infection, using the given strain distribution.
        '''
        is_a_sus = random.random() < self.sus/(self.sus + len(self.indvs))
        if is_a_sus:
            new_indv = individual(**self.indv_params)
            new_indv.setToMix(mix)
            self.indvs += [new_indv]
            for gnt in new_indv.getGenotypes(): self.inf[gnt] += 1
        else:
            random.shuffle(self.indvs)
            indv_to_infect = random.choice(self.indvs)
            init_gts = indv_to_infect.getGenotypes()
            indv_to_infect.infectSelfMult(mix)
            fin_gts = indv_to_infect.getGenotypes()
            if fin_gts != init_gts:
                gts_rmv = list(set(init_gts) - set(fin_gts))
                gts_add = list(set(fin_gts) - set(init_gts))
                for gt in gts_rmv: self.inf[gt] -= 1
                for gt in gts_add: self.inf[gt] += 1

    def updateSelBiases(self):
        if not self.do_sel_bias: return
        if self.is_vector: self.gnt_sel_bias = {sn: stats.mean(self.sel_bias_lst[sn]) for sn in self.sel_bias_lst}
        gsb_pre_parse = {sn: self.gnt_sel_bias[sn]/sum(self.gnt_sel_bias.values()) for sn in self.gnt_sel_bias}
        self.gnt_sel_bias: dict[str, float] = {}
        for sn in gsb_pre_parse:
            if sn == sn.upper(): self.gnt_sel_bias[sn] = gsb_pre_parse[sn]/(gsb_pre_parse[sn] + gsb_pre_parse[sn.lower()])
        for ind in self.indvs: ind.gnt_sel_bias = self.gnt_sel_bias.copy()
        self.indv_params['gnt_sel_bias'] = self.gnt_sel_bias

    def refresh(self):
        '''
        Filters out individuals that have been 'marked for death.'
        '''
        self.indvs = [ind for ind in self.indvs if not ind.marked_for_death]
        return self.indvs
        
    def printDatStr(self):
        '''
        Returns the population's information as a string.
        '''
        return ''.join([f'Population {self.pn} (total {self.tot_pop})\n',
                        f'S:\t{float2SN(self.sus, p=3)}\n',
                        f'\tI\n',
                        '\n'.join([f'{k}:\t{float2SN(self.inf[k])}' for k in self.inf]),
                        f'\n\tR\n',
                        '\n'.join([f'{k}:\t{float2SN(self.rec[k])}' for k in self.rec])])
    
    def printDat(self):
        '''
        Prints the population's information to the console.
        '''
        print(self.printDatStr())

    @property
    def is_dip(self):
        return not self.is_hap
    
    @is_dip.setter
    def is_dip(self, value: bool):
        self.is_hap = not value
    
    @property
    def individuals(self) -> list[individual]:
        return self.indvs
    
    @individuals.setter
    def individuals(self, value: list[individual]):
        self.indvs = value
    
    @property
    def tot_pop(self):
        I = 0
        if self.do_indvs: I = len(self.individuals)
        else: I = sum(self.inf.values())
        return int(self.sus + I + sum(self.rec.values()))

    def __str__(self):
        return self.pn
    
    def __repr__(self):
        return self.__str__()