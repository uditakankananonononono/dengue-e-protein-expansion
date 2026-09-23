#!/usr/bin/env python3
"""Build curated antibody-epitope residue sets for DENV-3/DENV-4 E protein.
Sources: Fibriansah 2015 Science (5J7 Table 1, PMC4346626); de Alwis 2012 PNAS (5J7 escape);
Sukupolvi-Petty 2013 JVI (DENV-4 DV4-E75 escape, PMC3754038); JVI 2007 5H2 escape K174;
PDB 4UT6 (EDE2 B7) and 4UTA (EDE1 C8) contact analysis (<5A, this repo code);
FLE canonical fusion loop 98-110. DENV-2 numbering mapped to DENV-3/4 refs by pairwise alignment."""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import json, re

def ref_e(dataset_dir):
    ref = SeqIO.read(f"{dataset_dir}/reference.fasta", "fasta")
    for line in open(f"{dataset_dir}/genome_annotation.gff3"):
        p = line.split("\t")
        if len(p) > 8 and p[2] == "gene" and p[8].strip() in ("gene_name=E", "E"):
            s, e = int(p[3]), int(p[4])
            return str(Seq(str(ref.seq[s-1:e])).translate())
    raise RuntimeError("E not found")

d2 = SeqIO.read("/tmp/denv2_ref.gb", "genbank")
d2e = str(Seq(str(d2.seq[935:2421])).translate())
d3e = ref_e("/tmp/nc_denv3"); d4e = ref_e("/tmp/nc_denv4")
print("E lengths: d2", len(d2e), "d3", len(d3e), "d4", len(d4e))

al = PairwiseAligner(); al.mode = "global"
al.match_score = 2; al.mismatch_score = -1; al.open_gap_score = -5; al.extend_gap_score = -0.5
def posmap(a, b):
    aln = al.align(a, b)[0]
    (a_blocks, b_blocks) = aln.aligned
    m = {}
    for (a0, a1), (b0, b1) in zip(a_blocks, b_blocks):
        for i in range(a1 - a0):
            m[int(a0) + i + 1] = int(b0) + i + 1
    return m
m23 = posmap(d2e, d3e); m24 = posmap(d2e, d4e)
print("map coverage d2->d3:", len(m23), "d2->d4:", len(m24))

EDE2_B7 = [68,69,70,71,72,73,74,82,97,98,99,101,102,103,104,105,113,152,153,154,155,156,157,245,246,247,248,249]
EDE1_C8 = [68,69,70,71,72,73,74,77,83,84,97,98,99,100,101,102,103,104,105,106,113,115,148,158,246,247,248,249,274,275,309,310,311,323,362]
FLE = list(range(98,111))
J7 = [50,51,52,53,54,55,58,73,74,101,106,123,126,128,130,131,133,134,148,196,198,200,201,223,224,227,274,276,307,308,309]
def mp(lst, m): return sorted({m[p] for p in lst if p in m})

epitopes = {
  "denv3": {
    "5J7_quaternary": {"residues": J7, "source": "Fibriansah 2015 Science PMC4346626 Table1 (DENV-3 native numbering)", "serotype_specific": True},
    "5J7_escape_site": {"residues": [269,270], "source": "de Alwis 2012 PNAS 1200566109 (Lys insertion Q269-N270 escape)", "serotype_specific": True},
    "EDE1_C8": {"residues": mp(EDE1_C8, m23), "source": "PDB 4UTA contact analysis <5A (DENV-2 numbering mapped)", "serotype_specific": False},
    "EDE2_B7": {"residues": mp(EDE2_B7, m23), "source": "PDB 4UT6 contact analysis <5A (DENV-2 numbering mapped)", "serotype_specific": False},
    "FLE_fusion_loop": {"residues": mp(FLE, m23), "source": "canonical fusion loop 98-110", "serotype_specific": False},
    "glycan_sites": {"residues": mp([67,153], m23), "source": "N67/N153 E glycosylation (EDE dependence)", "serotype_specific": False}
  },
  "denv4": {
    "escape_5H2": {"residues": [174], "source": "JVI 2007 5H2 escape K174 (DENV-4 native)", "serotype_specific": True},
    "escape_DV4-E75": {"residues": [330,361,364], "source": "Sukupolvi-Petty 2013 JVI PMC3754038 (G330E; T361A+V364I)", "serotype_specific": True},
    "EDE1_C8": {"residues": mp(EDE1_C8, m24), "source": "PDB 4UTA contact analysis <5A (DENV-2 numbering mapped)", "serotype_specific": False},
    "EDE2_B7": {"residues": mp(EDE2_B7, m24), "source": "PDB 4UT6 contact analysis <5A (DENV-2 numbering mapped)", "serotype_specific": False},
    "FLE_fusion_loop": {"residues": mp(FLE, m24), "source": "canonical fusion loop 98-110", "serotype_specific": False},
    "glycan_sites": {"residues": mp([67,153], m24), "source": "N67/N153 E glycosylation (EDE dependence)", "serotype_specific": False}
  }
}
json.dump(epitopes, open("/home/sandbox/work/denv3-4/results/epitopes.json","w"), indent=1)
for s in epitopes:
    for k,v in epitopes[s].items(): print(s, k, v["residues"])
