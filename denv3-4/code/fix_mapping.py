#!/usr/bin/env python3
"""Fix epitope mapping: align DENV-2 ref E to per-serotype E consensus from the
nextclade-aligned E translations (the exact coordinate frame of the variability analysis)."""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from collections import Counter
import json

def consensus(fasta):
    recs = list(SeqIO.parse(fasta, "fasta"))
    L = len(recs[0].seq)
    cons = []
    for i in range(L):
        c = Counter(str(r.seq[i]) for r in recs)
        if '-' in c: del c['-']
        if 'X' in c: del c['X']
        cons.append(c.most_common(1)[0][0] if c else 'X')
    return "".join(cons)

d2 = SeqIO.read("/tmp/denv2_ref.gb", "genbank")
d2e = open("/tmp/d2e.txt").read().strip()
cons = {"denv3": consensus("/tmp/ncout3/nextclade.cds_translation.E.fasta"),
        "denv4": consensus("/tmp/ncout4/nextclade.cds_translation.E.fasta")}
for s, c in cons.items():
    print(s, "consensus len", len(c), "N-term:", c[:25])

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
# anchor sanity checks
for s, m in maps.items():
    c = cons[s]
    checks = {p: (m.get(p), c[m[p]-1] if p in m else None) for p in [67, 101, 153]}
    print(s, "anchors d2pos->(mapped, aa):", checks, "coverage:", len(m))
json.dump({s: {str(k): v for k, v in m.items()} for s, m in maps.items()},
          open("/tmp/e_posmaps.json", "w"))
