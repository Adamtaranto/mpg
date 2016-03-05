def ntnum(nt):
    '''Nucleotide to number: A: 0, C: 1, G:2, T: 3.

    Case insensitive.
    '''
    n = (ord(nt) & 6) >> 1
    n ^= n>>1
    return n


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
