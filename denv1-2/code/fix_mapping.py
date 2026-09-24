#!/usr/bin/env python3
"""Map DENV-2 (NC_001474) E numbering onto the per-serotype consensus coordinate frames
(nextclade-aligned E translations) for DENV-1 and DENV-2."""
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from collections import Counter
import json

def consensus(fasta):
    recs = list(SeqIO.parse(fasta, "fasta"))
    L = max(len(r.seq) for r in recs)
    out = []
    for i in range(L):
        c = Counter(str(r.seq[i]) for r in recs if len(r.seq) > i)
        c.pop('-', None); c.pop('X', None)
        out.append(c.most_common(1)[0][0] if c else 'X')
    return "".join(out)

d2e = open("/tmp/d2e.txt").read().strip()
cons = {"denv1": consensus("/tmp/ncout1/nextclade.cds_translation.E.fasta"),
        "denv2": consensus("/tmp/ncout2/nextclade.cds_translation.E.fasta")}
al = PairwiseAligner(); al.mode = "global"
al.match_score = 2; al.mismatch_score = -1; al.open_gap_score = -8; al.extend_gap_score = -0.5
def posmap(a, b):
    aln = al.align(a, b)[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for i in range(int(a1) - int(a0)):
            m[int(a0) + i + 1] = int(b0) + i + 1
    return m
maps = {s: posmap(d2e, c) for s, c in cons.items()}
for s, m in maps.items():
    c = cons[s]
    print(s, "consensus len", len(c), "anchors 67/101/153:",
          {p: (m.get(p), c[m[p]-1] if p in m else None) for p in [67, 101, 153]}, "coverage:", len(m))
json.dump({s: {str(k): v for k, v in m.items()} for s, m in maps.items()},
          open("/tmp/e_posmaps_d12.json", "w"))
