cimport cython
ctypedef unsigned long long int u64

@cython.boundscheck(False)
cdef char numnt(u64 num):
    '''Number to nucleotide: A: 0, C: 1, G:2, T: 3.
    '''
    #if num == 0:
    #    return 'A'
    #elif num == 1:
    #    return 'C'
    #elif num == 2:
    #    return 'G'
    #else:
    #    return 'T'
    cdef const char *s = 'ACGT'
    return s[num]

def hash2kmer(h, k):
    kmer = []
    for x in range(k):
        kmer.append(numnt((h >> (2*x)) & 0x03))
    return ''.join(reversed(kmer))


def iter_kmers(seq, k):
    '''Iterator over hashed k-mers in a string DNA sequence.
    '''
    cdef u64 n
    cdef u64 bitmask = 2**(2*k)-1  # Set lowest 2*k bits
    cdef u64 h = 0

    # Pre-load the first k-1 nucleotides into the hash value
    # For each kmer's end nucleotide, bit-shift, add the end and yield
    for end in range(len(seq)):
        n = (ord(seq[end]) & 6) >> 1
        n ^= n>>1
        h = ((h << 2) | n) & bitmask
        if end >= k:
            yield h
