import itertools as itl

import numpy as np

from .counter import TransitionCounter
from .generator import MarkovGenerator
from .util import (
    iter_kmers,
    hash2kmer,
)

def test_iter_kmers():
    dbs = "AACAGATCCGCTGGTTA"
    k = 2
    counts = np.zeros(4**k)
    for kmer in iter_kmers(dbs, k):
        counts[kmer] += 1
    assert counts.sum() == len(dbs) - k + 1, counts.sum()
    assert (counts == 1).all(), counts

def test_hash2kmer():
    k = 2
    hashes = range(4**k)
    kmers = map(''.join, list(itl.product(list('ACGT'), repeat=k)))
    for hsh, mer in zip(hashes, kmers):
        h2k = hash2kmer(hsh, k)
        assert h2k == mer, (hsh, mer, h2k)

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


