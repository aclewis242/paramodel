'''
The file for the stochastic SIRS model class.
'''

from func_lib import *
from population import *
from allele import *
from random import shuffle
import pickle

class SIR:
	'''
	The stochastic SIRS model class.
	'''
	pn = ''			# Short population name
	pn_full = ''	# Full population name
	sn = 'init'		# Strain name (i.e. genotype)
	alleles = []	# List of allele objects
	bd = -1.0		# Birth/death rate
	ir = -1.0		# Infection rate (obsolete)
	rr = -1.0		# Recovery rate
	wi = -1.0		# Waning immunity rate
	itr: dict[population, float] = {}	# Dict of population objects to contact rates (usually only one entry)
	bds = 0.0		# Birth/death rate (for susceptibles specifically; value assigned in setRs)
	Rs = []			# Transition rate matrix
	Es = []			# List of events (3-element arrays that change S, I, R respectively)
	num_Es = -1		# Length of Es
	is_vector = False				# Whether or not this model describes a vector
	other_pop: population = None	# The other population this model interacts with
	
	def __init__(self, p0: population, **kwargs):
		'''
		Initialises the model with the given parameters.

		### Parameters
		- `p0`: Initial population, as a 3-element list (S, I, R)
		- `pn`: The name (short) of the population this model belongs to
		- `pn_full`: The full name of the population this model belongs to
		- `sn`: The name of the strain this model belongs to (i.e. genotype, in most cases)
		- `bd`: Birth/death rate (only relevant to vectors)
		- `ir`: Infection rate (obsolete)
		- `rr`: Recovery rate (only relevant to hosts)
		- `wi`: Waning immunity rate (only relevant to hosts)
		- `itr`: Interspecific transmission rates from this model to other populations (dict pop:float)
		- `is_vector`: Whether or not this model describes a vector population (bool)
		'''
		self.pop = p0
		self.__dict__.update(kwargs)
		self.pop.is_vector = self.is_vector
		self.pop.pn = self.pn
		self.pop.addStrain(self.sn)
		E1 = [1, 0, 0]	# Birth (used 3 times, one for each pop type)
		E2 = [-1, 1, 0] # Infection
		E3 = [0, -1, 1] # Recovery
		E4 = [-1, 0, 0] # Death of susceptible
		E5 = [0, -1, 0] # Death of infected
		E6 = [0, 0, -1] # Death of recovered
		E7 = [1, 0, -1] # Waning immunity
		self.Es = [E1, E2, E3, E4, E5, E6, E7, E1, E1, E1] # first E1, E2 are deprecated
		self.num_Es = len(self.Es)

	def setRs(self):
		'''
		Generates the different transition rates based on the model's parameters and population.
		'''
		# consider making this return the result instead of store it
		S = self.pop.sus
		R = self.pop.rec[self.sn]
		is_hap = self.pop.is_hap
		p2 = self.other_pop
		I_UW = 0
		I_WS = 0
		for ind in self.pop.individuals:
			if ind.genotype_freqs[self.sn]:
				I_WS += ind.correction_det(sn=self.sn)	# weighted according to how prevalent the strain is inside the indv
				if is_hap: continue
				I_UW += ind.correction_det()			# unweighted infections (simple 'yes/no' on strain presence)
		self.Rs = [ 0,				# deprecated (formerly births)
					0,				# deprecated (formerly intra-population infections)
					self.rr*I_WS,	# recoveries
					self.bds*S,		# susceptible deaths
					self.bd*I_UW,	# infected deaths
					self.bd*R,		# recovered deaths
					self.wi*R,		# waning immunity
					self.bds*S,		# susceptible births
					self.bd*I_UW,	# infected births (as in births from infecteds, not newly-born infecteds)
					self.bd*R,		# recovered births (again as in births from recovereds)
					self.itr[p2]*I_WS*(p2.sus+p2.getSusInfNum(self.sn))/p2.tot_pop]
									# interspecific contacts (cross-pop infections)

	def trans(self, idx: int, rpt: int=1):
		'''
		Effects the changes in the population dictated by the simulation.

		### Parameters
		- `idx`: The index of the desired event, corresponding with the order of `Rs`.
		- `rpt`: The number of times to repeat said event.
		'''
		pop = self.pop
		self_sn = self.sn
		pc_2_trans = self.pop.pc_to_transmit
		if idx >= self.num_Es:
			pop = self.other_pop
			idx = 1
			num_inf = 0
			num_mixes = 0
			num_loops = 0
			num_failed_trans = 0
			max_loops = 10000
			indvs_lst = self.pop.individuals
			while num_inf < rpt: # loops until the required number of events have occurred, or until the hard limit has been reached (unusual)
				for indv in indvs_lst:
					if not indv.genotype_freqs[self_sn]: continue
					elif indv.correction(sn=self_sn):
						num_inf += 1
						if indv.is_hap: # weird conditional structure is to save on time
							pop.infectMix(indv.infectMix(pc_2_trans))
							num_mixes += 1
						elif indv.is_mixed_vec:
							if indv.is_mixed:
								pop.infectMix(indv.infectMix(pc_2_trans))
								num_mixes += 1
							elif not indv.doesContactTransmit(): num_failed_trans += 1
						elif not indv.doesContactTransmit(): num_failed_trans += 1
					if num_inf >= rpt: break
				num_loops += 1
				if num_loops >= max_loops: break
				shuffle(indvs_lst) # slow but necessary to ensure fairness
			rpt -= (num_mixes + num_failed_trans)
		if rpt: pop.addPop(list(np.multiply(self.Es[idx], rpt)), self_sn, pc_2_trans)
	
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
		
		*Now obsolete. Largely a holdover from when alleles had macroscopic (population-level), not microscopic (individual-level),
		effects on the simulation.*

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
		
		*Now obsolete. Largely a holdover from when alleles had macroscopic (population-level), not microscopic (individual-level),
		effects on the simulation.*

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
	
	def save(self, s_dir: str=''):
		s_fn = f'{s_dir}/{self.fn}'
		os.makedirs(s_dir, exist_ok=True)
		f = open(s_fn, 'wb')
		pickle.dump(self, f)
		f.close()

	def load(path: str) -> 'SIR':
		if not os.path.exists(path): print(f'failed to load model from path {path}'); return
		f = open(path, 'rb')
		mdl = pickle.load(f)
		f.close()
		return mdl
	
	@property
	def fn(self): # strange naming scheme is due to pickle pathing being case-insensitive (alleles don't work)
		disc = 'et'
		if self.sn == self.sn.upper(): disc = 'm'
		elif self.sn == self.sn.lower(): disc = 'w'
		if len(self.sn) == 2: disc = f'h{disc}'
		return f'{self.pn}_{disc}.sir'

	def __str__(self):
		return f'population {self.pn}, strain {self.sn}'
	
	def __repr__(self):
		return self.__str__()