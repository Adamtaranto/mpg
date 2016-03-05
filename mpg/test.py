
def test_ntnum():
    truth = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    for a in 'ACGTacgt':
        assert ntnum(a) == truth.get(a.lower(), 0), (a, ntnum(a))


def test_iter_kmers():
    dbs = "AACAGATCCGCTGGTTA"
    k = 2
    counts = np.zeros(4**k)
    for kmer in iter_kmers(dbs, k):
        counts[kmer] += 1
    assert counts.sum() == len(dbs) - k + 1, counts.sum()
    assert (counts == 1).all(), counts

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


