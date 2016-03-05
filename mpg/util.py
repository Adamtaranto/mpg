from ._util import (
    hash_to_kmer,
    iter_kmers,
)


def seq2fa(name, seq, linelen=80):
    lines = ['>{}'.format(name), ]

    for start in range(0, len(seq), linelen):
        lines.append(seq[start:start + linelen])

    return '\n'.join(lines) + '\n'
