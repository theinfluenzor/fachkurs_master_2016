import random
import sys

class BioMolecule(object):
    """
    A generic molecule that has basic attributes like id, name and
    mass.
    """
    def __init__(self, id, name, mass=0):
        self._id = id
        self.name = name
        self.mass = mass

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        self._name = value

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        # if self.mass and \
        #    not isinstance(value, float) or \
        #    not isinstance(value, int):
        #     raise Exception("mass must be numeric")
#        else:
        self._mass = value


class Polymer(BioMolecule):
    """
    A polymer molecule that has a sequence attribute which is
    accessible via indexing the object.
    """
    def __init__(self, id, name, sequence, mass=0):
        super(Polymer, self).__init__(id, name, mass)
        self._sequence = sequence

    def __getitem__(self, value):
        return self.sequence[value]

    def __setitem__(self, key, value):
        self.sequence[key] = value

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise Exception("sequence must be a string")
            # TODO: check for valid nucleotides here
        self._sequence = value.upper()
        self.calculate_mass()


class MRNA(Polymer):
    def __init__(self, id, name, sequence, mass=0):
        super(MRNA, self).__init__(id, name, sequence, mass)
        self.binding = [0]*(len(sequence)/3)

    def calculate_mass(self):
        self.mass = 0
        NA_mass = {'A': 1.0, 'U': 2.2, 'G':2.1, 'C':1.3}
        for na in self.sequence:
            self.mass += NA_mass[na]


class Protein(Polymer):
    number_of_proteins = 0

    def __init__(self, id, name, sequence, mass=0):
        super(Protein, self).__init__(id, name, sequence, mass)
        self.number_of_proteins += 1

    def __add__(self, AS):
        self.sequence += AS 

    def calculate_mass(self):
        AA_mass = {"A": 89.0929, "R": 175.208, "N": 132.118, "D": 132.094, "C": 121.158, "Q": 146.144,
                    "E": 146.121, "G": 75.0664, "H":155.154, "I":131.172, "L": 131.172, "K": 147.195,
                    "M": 149.211, "F": 165.189, "P": 115.13, "S": 105.092, "T": 119.119, "W": 204.225,
                    "Y":181.188, "V":117.146}
        for aa in self.sequence:
            self.mass += AA_mass[aa]
   

class Ribosome(BioMolecule):
    """
    """
    code = dict([('UCA','S'), ('UCG','S'), ('UCC','S'), ('UCU','S'),
                 ('UUU','F'), ('UUC','F'), ('UUA','L'), ('UUG','L'),
                 ('UAU','Y'), ('UAC','Y'), ('UAA','*'), ('UAG','*'),
                 ('UGU','C'), ('UGC','C'), ('UGA','*'), ('UGG','W'),
                 ('CUA','L'), ('CUG','L'), ('CUC','L'), ('CUU','L'),
                 ('CCA','P'), ('CCG','P'), ('CCC','P'), ('CCU','P'),
                 ('CAU','H'), ('CAC','H'), ('CAA','Q'), ('CAG','Q'),
                 ('CGA','R'), ('CGG','R'), ('CGC','R'), ('CGU','R'),
                 ('AUU','I'), ('AUC','I'), ('AUA','I'), ('AUG','M'),
                 ('ACA','T'), ('ACG','T'), ('ACC','T'), ('ACU','T'),
                 ('AAU','N'), ('AAC','N'), ('AAA','K'), ('AAG','K'),
                 ('AGU','S'), ('AGC','S'), ('AGA','R'), ('AGG','R'),
                 ('GUA','V'), ('GUG','V'), ('GUC','V'), ('GUU','V'),
                 ('GCA','A'), ('GCG','A'), ('GCC','A'), ('GCU','A'),
                 ('GAU','D'), ('GAC','D'), ('GAA','E'), ('GAG','E'),
                 ('GGA','G'), ('GGG','G'), ('GGC','G'), ('GGU','G')])

    def __init__(self, id, name):
        super(Ribosome, self).__init__(id, name)
        self.bound = False
        self.position = None

    def initiate(self, mrna):
        if not self.bound and not mrna.binding[0]:  #  no mrna bound yet and target mrna still free at pos 0
            self.bound = mrna
            self.nascent_prot = Protein(self.bound.id,
                                        "Protein_{0}".format(self.bound.id),
                                        "")
            self.position = 0
            self.bound.binding[0] = 1
            
    def elongate(self):
        if not self.bound: # can't elongate
            return False
        codon = self.bound[self.position:self.position+3]
        aa = self.code[codon]

        if aa == "*": # terminate at stop codon
            return self.terminate()
            
        if not self.bound.binding[self.position/3 + 1]: # if the next rna position is free
            self.bound.binding[self.position/3] = 0
            self.bound.binding[self.position/3+1] = 1
            self.position += 3
            self.nascent_prot + aa
        return 0

    def terminate(self):
        self.bound.binding[self.position/3] = 0 # bound mRNA
        self.bound = False
        return self.nascent_prot
        

class Cell(object):
    def __init__(self):
        self.ribosomes = [Ribosome(i, 'Ribo_{0}'.format(i)) for i in xrange(200)]
        self.mrnas = [MRNA(i, 'MRNA_{0}'.format(i), "UUUUUUUUUUAA") for i in xrange(20)]
        self.proteins = [[] for x in xrange(20)]

    def step(self):
        for r in self.ribosomes:
            if not r.bound:
                r.initiate(self.mrnas[random.randint(0,len(self.mrnas)-1)])
            else:
                prot = r.elongate()
                if prot:
                    self.proteins[prot.id].append(prot)

    def simulate(self, steps, p=True):
        for s in xrange(steps):
            self.step()
            if p:
                print [len(x) for x in self.proteins]
            
if __name__ == "__main__":
    c = Cell()
    c.simulate(10)
