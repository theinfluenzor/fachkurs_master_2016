import processes as proc
import molecules as mol

class Cell(object):
    def __init__(self):
        self.states = {}
        
        self.ribosomes = {'Ribo_{0}'.format(i): mol.Ribosome(i, 'Ribo_{0}'.format(i)) for i in xrange(10)}
        self.mrnas = {'MRNA_{0}'.format(i): mol.MRNA(i, 'MRNA_{0}'.format(i), "UUUUUUUUUUAA") for i in xrange(20)}
        self.proteins = [[] for x in xrange(20)]
        self.states.update(self.ribosomes)
        self.states.update(self.mrnas)
        
        translation = proc.Translation(1, "Translation")
        translation.set_states(self.mrnas.keys(), self.ribosomes.keys())
        self.processes = {"Translation":translation}

    def step(self):
        for p in self.processes:
            self.processes[p].update(self)

    def simulate(self, steps, log=True):
        for s in xrange(steps):
            self.step()
            if log:
                # print the count of each protein to the screen
                print '\r{}'.format([len(self.states[x]) for x in self.states.keys() if "Protein_" in x]),
            
if __name__ == "__main__":
    c = Cell()
    c.simulate(1000, log=True)
