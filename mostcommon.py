from collections import Counter
import screed
import sys

K=3
c = Counter()
with screed.open(sys.argv[1]) as fh:
    for rec in fh:
        s = rec.sequence
        l = len(s)
        for start in range(l - K + 1):
            c[s[start:start+K]] += 1

for kmer, cnt in c.most_common(100):
    print kmer, cnt
