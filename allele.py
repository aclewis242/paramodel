'''
The file containing the class for alleles.
'''

class allele:
	'''
	The class for alleles. Describes the phenotype observed from its expression and contains the character representing it in text.
	'''
	char = ''
	fav_pop = ''
	unf_pop = ''
	param = ''
	fac = 0.0
	sel_advs: dict[str, float] = {} 			# specifically for mutated allele (capital)
	transm_probs: dict[str, float] = {} 		# also specifically for mutated allele
	base_transm_probs: dict[str, float] = {} 	# specifically for wild allele (lowercase)

	def __init__(self, **kwargs):
		'''
		Initialises the allele.

		### Parameters
		- `char`: The character representing the allele in text
		- `fav_pop`: The name of the population it's well-adapted to *(obsolete)*
		- `unf_pop`: The name of the population it's poorly-adapted to *(obsolete)*
		- `param`: The name of the trait it affects *(obsolete)*
		- `fac`: The numerical factor by which the trait is affected. 'Good' traits should be positive, 'bad' traits should be negative
			(its value should be between -1 and 1, not inclusive) *(obsolete)*
		'''
		self.__dict__.update(kwargs)
		self.char = self.char.upper()
		self.sel_advs: dict[str, float] = {}
		self.transm_probs: dict[str, float] = {}
		self.base_transm_probs: dict[str, float] = {}
	
	@property
	def locus(self):
		'''
		The 'ID' of the locus this allele is at. (This is the lowercase of the allele's character.)
		'''
		return self.char.lower()
	
	def __str__(self):
		return self.char
	
	def __repr__(self):
		return self.__str__()