from collections import deque

import screed
import numpy as np
from numpy import random
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg
import yaml


def seq2fa(name, seq, linelen=80):
    lines = ['>{}'.format(name),]

    for start in range(0, len(seq), linelen):
        lines.append(seq[start:start+linelen])

    return '\n'.join(lines) + '\n'


def ntnum(nt):
    '''Nucleotide to number: A: 0, C: 1, G:2, T: 3.

    Case insensitive.
    '''
    n = (ord(nt) & 6) >> 1
    n ^= n>>1
    return n


def test_ntnum():
    truth = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    for a in 'ACGTacgt':
        assert ntnum(a) == truth.get(a.lower(), 0), (a, ntnum(a))


def numnt(num):
    '''Number to nucleotide: A: 0, C: 1, G:2, T: 3.
    '''
    return "ACGT"[num]

def hash2kmer(h, k):
    kmer = []
    for x in range(k):
        kmer.append(numnt((h >> (2*x)) & 0x03))
    return ''.join(reversed(kmer))


def iter_kmers(seq, k):
    '''Iterator over hashed k-mers in a string DNA sequence.
    '''
    bitmask = 2**(2*k)-1  # Set lowest 2*k bits
    h = 0

    # Pre-load the first k-1 nucleotides into the hash value
    for nt in seq[:k-1]:
        h = (h << 2) | ntnum(nt)
    # For each kmer's end nucleotide, bit-shift, add the end and yield
    for end in range(k-1, len(seq)):
        h = ((h << 2) | ntnum(seq[end])) & bitmask
        yield h


def test_iter_kmers():
    dbs = "AACAGATCCGCTGGTTA"
    k = 2
    counts = np.zeros(4**k)
    for kmer in iter_kmers(dbs, k):
        counts[kmer] += 1
    assert counts.sum() == len(dbs) - k + 1, counts.sum()
    assert (counts == 1).all(), counts


class TransitionCounter(object):

    def __init__(self, k, alphabet=set("ACGT")):
        self.bitmask = 2**(2*k)-1  # Set lowest 2*k bits
        self.k = k
        self.alphabet = set(alphabet)
        self.n = len(alphabet) ** k
        self.transition_counts = np.zeros((self.n, len(alphabet)))
        self._transitions = None
        self._P = None

    def _clear(self):
        self._transitions = None
        self._P = None

    def save(self, filename):
        data = {
            'alphabet': list(sorted(self.alphabet)),
            'k': self.k,
            'transitions': self.transitions.tolist(),
        }
        with open(filename, 'w') as fh:
            yaml.dump(data, fh)

    def load(self, filename):
        with open(filename, 'r') as fh:
            data = yaml.load(fh)
        datakeys = list(sorted(data.keys()))
        if datakeys != ['alphabet', 'k', 'transitions']:
            exc = ValueError("Data file is not valid")
            exc.keys = datakeys
            raise exc
        self.__init__(data['k'], data['alphabet'])
        self._transitions = np.array(data['transitions'])
        self.transition_counts = np.array(data['transitions'])

    def consume(self, sequence):
        self._clear()
        for subseq in sequence.split('N'):
            if not subseq or len(subseq) < self.k:
                continue
            kmers = iter_kmers(subseq, self.k)
            # pop the first kmer
            fr = next(kmers)
            for to in kmers:
                self.transition_counts[fr, to & 3] += 1
                fr = to

    def consume_file(self, filename):
        self._clear()
        with screed.open(filename) as sequences:
            for seq in sequences:
                self.consume(seq['sequence'])

    @property
    def transitions(self):
        if self._transitions is not None:
            return self._transitions
        transitions = self.transition_counts
        transitions /= transitions.sum(1)[:, np.newaxis]
        self._transitions = transitions
        return transitions

    @property
    def P(self):
        if self._P is not None:
            return self._P
        sparse_P = sparse.lil_matrix((self.n, self.n))
        num_kmers, alpha_size = self.transitions.shape
        for fr in range(num_kmers):
            for a in range(alpha_size):
                to = (fr << 2 | a) & self.bitmask
                sparse_P[fr, to] = self.transitions[fr, a]
        self._P = sparse_P
        return sparse_P

    @property
    def steady_state(self):
        v, w = linalg.eigs(self.P.transpose(), which='LR')
        ssf = np.real(w[:, v.argmax()])
        ssf /= ssf.sum()
        return ssf


def test_transition_counter_consume():
    dbs = 'AAACAAGAATACCACGACTAGCAGGAGTATCATGATTCCCGCCTCGGCGTCTGCTTGGGTGTTTAA'
    t = TransitionCounter(2)
    t.consume(dbs)
    counts = t.transition_counts
    assert (counts == 1).all(), counts
    P = t.transitions
    assert (P.sum(1) == 1).all(), P.sum(1)
    pi = t.steady_state
    assert pi.sum() == 1, pi.sum()
    assert np.allclose(pi.dot(t.P.toarray()),  pi), pi


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
        sequence = list(hash2kmer(prev_mer, self.k))

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

CLI = '''
USAGE:
    mpg [options] <reference> ...

OPTIONS:
    -l LENGTH       Length of sequence to simulate [default: 0].
    -k ORDER        Markovian Order [default: 1].
    -s SEED         Random seed for generation. (uses /dev/urandom if unset)
    -d DUMPFILE     Dump data to DUMPFILE.
    -r              Reference file is a Yaml dump file.
'''

if __name__ == '__main__':
    from docopt import docopt
    from sys import stdout, stderr
    opts = docopt(CLI)
    k = int(opts['-k'])
    l = int(opts['-l'])

    t = TransitionCounter(k)
    if opts['-r']:
        filename = opts['<reference>'][0]
        print('Loading reference dump from "{}"'.format(filename),
              file=stderr)
        t.load(filename)
    else:
        print('Inferring transition matrix from reference sequences',
              file=stderr)
        for fn in opts['<reference>']:
            t.consume_file(fn)
            print('Consumed', fn, file=stderr)

    if opts['-d']:
        filename = opts['-d']
        print('Saving reference dump to "{}"'.format(filename), file=stderr)
        t.save(filename)

    m = MarkovGenerator(t, seed=opts['-s'])
    print('Initialised Markov Generator', file=stderr)
    if l > 0:
        print('Generating sequence of {} bases'.format(l), file=stderr)
        print(seq2fa('genseq', m.generate_sequence(l)))
