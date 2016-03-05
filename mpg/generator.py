from numpy import random

from .util import (
    hash_to_kmer,
)

class MarkovGenerator(object):

    def __init__(self, transcount, seed=None):
        self.P = transcount.transitions
        self.pi = transcount.steady_state
        self.k = transcount.k
        self.bitmask = 2**(2*self.k)-1
        self.rand = random.RandomState()
        self.rand.seed(seed)

    def generate_sequence(self, length, seed=None):
        if seed is not None:
            self.rand.seed(seed)

        prev_mer = self.rand.choice(self.pi.size, p=self.pi)
        sequence = list(hash_to_kmer(prev_mer, self.k))

        for x in range(self.k, length):
            # Emission probs for the previous kmer
            p = self.P[prev_mer]
            nt = self.rand.choice(p.size, p=p)
            sequence.append(numnt(nt))
            # Add the nucleotide to the previous hash using bit ops
            # equiv. to x = x[1:] + nt
            prev_mer <<= 2
            prev_mer |= nt
            prev_mer &= self.bitmask

        return ''.join(sequence)
