import random
import molecules


class Process(object):
    """
    Parent for all cellular processes.
    """
    def __init__(self, id, name):
        self.id = id
        self.name = name

        self.enzyme_ids = []
        self.substrate_ids = []

    def set_states(self, substrate_ids, enzyme_ids):
        self.enzyme_ids = enzyme_ids
        self.substrate_ids = substrate_ids

    def update(self, cell):
        """
        Has to be implemented by child class.
        """
        pass


class Translation(Process):
    """
    Translation is instantiated in the Cell to produce proteins.

    Defines Translation process. It iterates over all ribosomes and decides what
    they should do. They either bind to mRNA or elongate/terminate a protein if
    they are already bound.

    """
    def __init__(self, id, name):
        super(Translation, self).__init__(id, name)

    def update(self, cell):
        for ribosome_id in self.enzyme_ids:
            prot = None
            ribosome = cell.states[ribosome_id]
            if not ribosome.bound:
                random_mrna = cell.states[self.substrate_ids[random.randint(0,len(self.substrate_ids)-1)]]
                ribosome.initiate(random_mrna)
            else:
                prot = ribosome.elongate
            if isinstance(prot, molecules.Protein):
                if prot.id in cell.states:
                    cell.states[prot.id].append(prot)
                else:
                    cell.states[prot.id] = [prot]