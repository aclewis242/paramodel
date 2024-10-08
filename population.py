from func_lib import *
from individual import *
from random import shuffle, choice

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
    indv_params = {}
    rng: np.random.Generator = None
    pc_to_transmit = 0
    all_sel_bias: dict[str, float] = {}
    all_transm_probs: dict[str, float] = {}
    init_pop = 0
    times: list[float] = []
    gnts: list[str] = []
    num_gnts = 0
    main_all_char: str = ''

    def __init__(self, p0: list[int], pn: str='', isn: str='init', **kwargs):
        '''
        Initialises the population.

        ### Parameters
        - `p0`: A 3-element integer list of the starting population amounts (S, I, R)
        - `pn`: The population's name
        - `isn`: The name of the initial strain
        - `is_vector`: Whether or not this population is a vector
        - `is_hap`: Whether or not this population is haploid
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
        self.indv_params['rng'] = self.rng
        self.indv_params['pc_flt'] = float(self.indv_params['pc'])
        self.all_sel_bias: dict[str, float] = {}
        self.all_transm_probs: dict[str, float] = {}
        self.times: list[float] = [0 for i in range(15)]
        self.gnts: list[str] = []
    
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
        if weight:
            infs = dict.fromkeys(infs, 0)
            for ind in self.indvs:
                for gt in ind.genotype_freqs: infs[gt] += ind.genotype_freqs[gt]/ind.pc_flt
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
        tm = time()
        sn = self.match(sn)
        # old_data = self.printDatStr()
        # old_indvs_len = len(self.indvs)
        # old_sus = self.sus
        # old_real_infs = dict.fromkeys(self.inf, 0)
        # for ind in self.indvs:
        #     for gt in ind.getGenotypes(): old_real_infs[gt] += 1
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
        self.times[0] += time() - tm
        tm = time()
        S = int(p[0])
        I = int(p[1])
        R = int(p[2])
        sus_new = self.sus
        inf_new = self.inf.copy()
        rec_new = self.rec.copy()
        self.times[1] += time() - tm

        tm = time()
        sus_new += S
        inf_new[sn] += I
        inf_indvs = []
        i_change = 0
        to_change = {}
        dead_gnts = []
        tba: list[individual] = []
        # mod_inds: list[individual] = []
        # pre_gnts: list[list[str]] = []
        self.times[3] += time() - tm
        if I > 0:
            tm = time()
            sus_indvs = self.getSusInf(sn)
            strain_indvs = {sn_i: self.getSusInf(sn_i, is_present=True, indvs_lst=sus_indvs) for sn_i in self.inf}
            sus_pop_nums = [self.sus] + [len(strain_indvs[sn_i]) if sn != sn_i else 0 for sn_i in self.inf]
            sus_pop_wgts = normalise(sus_pop_nums)
            self.times[2] += time() - tm
            tm = time()
            i_change, to_change = self.getChanges(I, sus_pop_wgts)
            tba = self.makeIndvs(sn, i_change)
            self.times[4] += time() - tm
            for strn in to_change:
                tm = time()
                sus_new += to_change[strn]
                self.times[5] += time() - tm
                tm = time()
                indvs_to_infect = strain_indvs[strn]
                self.times[6] += time() - tm
                tm = time()
                shuffle(indvs_to_infect)
                self.times[7] += time() - tm
                for ind in indvs_to_infect[:to_change[strn]]:
                    tm = time()
                    dead_gnts = ind.infectSelf(pc_trans, sn, do_return=True)
                    for d_g in dead_gnts: inf_new[d_g] -= 1
                    # mod_inds += [ind]
                    self.times[8] += time() - tm
            tm = time()
            self.indvs += tba
            self.times[9] += time() - tm
        elif I < 0:
            # self.refresh()
            tm = time()
            # random.shuffle(self.indvs)
            inf_indvs = self.getSusInf(sn, is_present=True)
            self.times[10] += time() - tm
            tm = time()
            shuffle(inf_indvs)
            self.times[11] += time() - tm
            for ind in inf_indvs[:-I]:
                # if not R:
                #     f = open('inf_events_raw.dat', 'a')
                #     f.write(f'death {ind.genotype_freqs}\n')
                #     f.close()
                tm = time()
                ind.marked_for_death = True
                other_gts = list(set(ind.getGenotypes()) - set([sn]))
                for gt in other_gts: inf_new[gt] -= 1
                self.times[12] += time() - tm
            tm = time()
            self.refresh()
            self.times[13] += time() - tm
        tm = time()
        rec_new[sn] += R

        self.sus = sus_new
        self.inf = inf_new
        self.rec = rec_new
        self.times[14] += time() - tm
        # self.refresh()
        # f = open('last_event.dat', 'w')
        # f.write(f'\n-----\nadded {p} with strain {sn}\nold data:\n')
        # f.write(old_data)
        # f.write('\nnew data:\n')
        # f.write(self.printDatStr())
        # f.close()

        # real_infs = dict.fromkeys(self.inf, 0)
        # for ind in self.indvs:
        #     for gt in ind.getGenotypes(): real_infs[gt] += 1
        # if (self.init_pop and not self.is_vector) and (self.tot_pop != self.init_pop or self.inf != real_infs):
        #     print(f'p: {p}, sn: {sn}')
        #     print('OLD:')
        #     print(old_data)
        #     print(f'len old indvs: {old_indvs_len}; old sus: {old_sus}')
        #     print('NEW:')
        #     self.printDat()
        #     print(f'len new indvs: {len(self.indvs)}; new sus: {self.sus}')
        #     print(sum(self.rec.values()))
        #     [print(ind.genotype_freqs) for ind in inf_indvs]
        #     print(f'I_corr: {I_corr}')
        #     print(f'i_change: {i_change}')
        #     print(f'to_change: {to_change}')
        #     print(f'old_real_infs: {old_real_infs}')
        #     print(f'real_infs: {real_infs}')
        #     print(f'self.inf: {self.inf}')
        #     print(f'dead_gnts: {dead_gnts}')
        #     print(f'tba gnts: {[ind.getGenotypes() for ind in tba]}')
        #     print(f'mod_ind gnts: {[ind.getGenotypes() for ind in mod_inds]}')
        #     print(f'pre_gnts: {pre_gnts}')
        #     print(f'sus_indvs gnts: {[ind.getGenotypes() for ind in sus_indvs]}')
        #     exit()
    
    def getChanges(self, pop_num: int, weights: list[float]) -> tuple[int, dict[str, int]]:
        '''
        Distributes infections between full-susceptibles and infected-susceptibles.
        '''
        to_change_lst = self.rng.multinomial(pop_num, weights)
        return to_change_lst[0], dictify(self.inf.keys(), to_change_lst[1:])
    
    def makeIndvs(self, sn: str, num_indvs: int):
        '''
        Make the given number of individuals with the given strain.
        '''
        return [individual(gnts=self.gnts, gnt=sn, **self.indv_params) for i in range(int(num_indvs))]

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
    
    def getSusInf(self, sn: str, is_present: bool=False, indvs_lst: list[individual]=None):
        '''
        Gets list of infected individuals either with strain `sn` (`is_present`) or susceptible to infection by the strain
        (`not is_present`). A list to draw from can be provided, if desired (otherwise, it'll draw from `self.indvs`).
        '''
        is_present -= 1
        if indvs_lst is None: indvs_lst = self.indvs
        return [ind for ind in indvs_lst if bool(ind.genotype_freqs[sn])+is_present]
    
    def getSusInfNum(self, sn: str, is_present: bool=False, indvs_lst: list[individual]=[]):
        '''
        Gets the number of infected individuals either with strain `sn` (`is_present`) or susceptible to infection by the strain
        (`not is_present`). A list to draw from can be provided, if desired (otherwise, it'll draw from `self.indvs`).
        '''
        rv = 0
        sn = self.match(sn)
        if not indvs_lst: indvs_lst = self.indvs
        is_present -= 1
        for ind in indvs_lst:
            if bool(ind.genotype_freqs[sn])+is_present: rv += 1
        return rv

    def infectMix(self, mix: dict[str, int]):
        '''
        Perform a mixed infection, using the given strain distribution.
        '''
        is_empty = True
        for gt in mix:
            if mix[gt]: is_empty = False
        if is_empty: return
        is_a_sus = random() < self.sus/(self.sus + len(self.indvs))
        if is_a_sus:
            new_indv = individual(gnts=self.gnts, **self.indv_params)
            new_indv.setToMix(mix)
            self.indvs += [new_indv]
            for gnt in new_indv.genotype_freqs:
                if new_indv.genotype_freqs[gnt]: self.inf[gnt] += 1
            self.sus -= 1
        else:
            # random.shuffle(self.indvs)
            indv_to_infect = choice(self.indvs)
            init_gts = indv_to_infect.getGenotypes()
            indv_to_infect.infectSelfMult(mix)
            fin_gts = indv_to_infect.getGenotypes()
            if fin_gts != init_gts: # can likely be made meaningfully faster
                gts_rmv = list(set(init_gts) - set(fin_gts))
                gts_add = list(set(fin_gts) - set(init_gts))
                for gt in gts_rmv: self.inf[gt] -= 1
                for gt in gts_add: self.inf[gt] += 1
    
    def updateSelBiases(self, alleles: list[allele]):
        '''
        Updates the population's selection (as well as transmission) biases based on the properties of the given alleles.
        '''
        self.all_sel_bias: dict[str, float] = {}
        self.all_transm_probs: dict[str, float] = {}
        for a in alleles:
            self.all_sel_bias[a.char] = a.sel_advs[self.pn]
            self.all_transm_probs[a.char] = a.transm_probs[self.pn]
            self.all_transm_probs[a.char.lower()] = a.base_transm_probs[self.pn]
            self.main_all_char = a.char.upper()
        for ind in self.indvs:
            ind.all_sel_bias = self.all_sel_bias.copy()
            ind.all_transm_probs = self.all_transm_probs.copy()
            ind.main_all_char = self.main_all_char
        self.indv_params['all_sel_bias'] = self.all_sel_bias
        self.indv_params['all_transm_probs'] = self.all_transm_probs
        self.indv_params['main_all_char'] = self.main_all_char

    def getGntDist(self):
        '''
        Gets all the individuals' (capital/mutated) allele frequencies.
        '''
        return [ind.getAlleleFreqs()/ind.num_genes for ind in self.indvs]

    def refresh(self, update: bool=True):
        '''
        Filters out individuals that have been 'marked for death,' and (optionally) updates the macroscopic infection tracking to make sure
        it aligns with the microscopic state of the simulation.
        '''
        self.indvs = [ind for ind in self.indvs if not ind.marked_for_death] # consider trying to fold the two for loops here together
        if update:
            self.inf = dict.fromkeys(self.inf, 0)
            for ind in self.indvs:
                gtfs = ind.genotype_freqs
                for gt in gtfs:
                    if gtfs[gt]: self.inf[gt] += 1
        
    def printDatStr(self):
        '''
        Returns the population's information as a string.
        '''
        return ''.join([f'Population {self.pn} (total {self.tot_pop})\n',
                        f'S:\t{float2SN(self.sus, p=3)}\n',
                        f'\tI\n',
                        '\n'.join([f'{k}:\t{float2SN(self.inf[k])}' for k in self.inf]),
                        f'\nTotal:\t{self.num_inf}',
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
        if self.is_hap and self.init_pop: return self.init_pop
        return int(self.sus + self.num_inf + sum(self.rec.values()))
    
    @property
    def num_inf(self):
        return len(self.individuals)

    def __str__(self):
        return self.pn
    
    def __repr__(self):
        return self.__str__()