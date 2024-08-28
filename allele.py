class allele:
    '''
    The class for alleles. Describes the phenotype observed from its expression and contains the character representing it on the genome.
    '''
    char = ''
    fav_pop = ''
    unf_pop = ''
    param = ''
    fac = 0.0

    def __init__(self, **kwargs):
        '''
        Initialises the allele.

        ### Parameters
        - `char`: The character representing the allele on the genome
        - `fav_pop`: The name of the population it's well-adapted to
        - `unf_pop`: The name of the population it's poorly-adapted to
        - `param`: The name of the trait it affects
        - `fac`: The numerical factor by which the trait is affected. 'Good' traits should be positive, 'bad' traits should be negative
            (its value should be between -1 and 1, not inclusive)
        '''
        self.__dict__.update(kwargs)
        self.char = self.char.upper()
    
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