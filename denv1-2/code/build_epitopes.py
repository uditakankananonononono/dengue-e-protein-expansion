#!/usr/bin/env python3
"""Curated antibody-epitope residue sets for DENV-1/DENV-2 E protein.
Sources: de Alwis 2012 PNAS 1200566109 (1F4 escape K47E/G274E, DENV-1 native; 2D22 escape
R323G + associated H282/D362, DENV-2 native); Shrestha 2010 JVI / Diamond lab (DENV1-E106
DIII epitope, PMC3312934 escape T329A); PDB 4C2I (1F4-DENV1 virion) and 4UIF (2D22-DENV2
virion) contact analysis <=6A (this repo, code/pdb_contacts.py); PDB 4UTA/4UT6 EDE contacts
(<5A, DENV-2 numbering; mapped to DENV-1 by pairwise alignment); canonical FLE 98-110."""
import json

maps = json.load(open("/tmp/e_posmaps_d12.json"))
EDE2_B7 = [68,69,70,71,72,73,74,82,97,98,99,101,102,103,104,105,113,152,153,154,155,156,157,245,246,247,248,249]
EDE1_C8 = [68,69,70,71,72,73,74,77,83,84,97,98,99,100,101,102,103,104,105,106,113,115,148,158,246,247,248,249,274,275,309,310,311,323,362]
FLE = list(range(98,111)); GLY = [67,153]
# PDB contact analysis outputs (code/pdb_contacts.py), <=6A all-atom:
F1F4 = [51,155,156,162,163,164,165,171,176,177]          # 4C2I, EDI part (hinge 274/275 kept as GT only)
F2D22 = [67,68,69,70,71,72,73,101,102,103,104,152,153]    # 4UIF fusion-loop/DII footprint
E106 = [310,325,328,330,361,362,364,385]                   # Shrestha 2010 (T329 escape excluded, kept as GT)
def mp(lst, s): return sorted({maps[s][str(p)] for p in lst if str(p) in maps[s]})

epitopes = {
  "denv1": {
    "1F4_structural": {"residues": F1F4, "source": "PDB 4C2I contact analysis <=6A (this work); Fibriansah 2014 Science 4C2I", "role": "feature"},
    "E106_DIII": {"residues": E106, "source": "Shrestha 2010 JVI PMC3312934 (K310,P325,G328,D330,K361,E362,P364,K385; T329 escape kept as ground truth)", "role": "feature"},
    "EDE1_C8": {"residues": mp(EDE1_C8,"denv1"), "source": "PDB 4UTA contact analysis <5A (DENV-2 numbering mapped)", "role": "feature"},
    "EDE2_B7": {"residues": mp(EDE2_B7,"denv1"), "source": "PDB 4UT6 contact analysis <5A (DENV-2 numbering mapped)", "role": "feature"},
    "FLE_fusion_loop": {"residues": mp(FLE,"denv1"), "source": "canonical fusion loop 98-110 (mapped)", "role": "feature"},
    "glycan_sites": {"residues": mp(GLY,"denv1"), "source": "N67/N153 E glycosylation (EDE dependence)", "role": "feature"},
    "escape_1F4": {"residues": [47,274], "source": "de Alwis 2012 PNAS 1200566109 (K47E, G274E; DENV-1 native)", "role": "ground_truth"},
    "escape_E106": {"residues": [329], "source": "PMC3312934 (T329A; DENV-1 native)", "role": "ground_truth"}
  },
  "denv2": {
    "2D22_structural": {"residues": F2D22, "source": "PDB 4UIF contact analysis <=6A (this work); de Alwis 2015 Science", "role": "feature"},
    "EDE1_C8": {"residues": mp(EDE1_C8,"denv2"), "source": "PDB 4UTA contact analysis <5A (DENV-2 numbering)", "role": "feature"},
    "EDE2_B7": {"residues": mp(EDE2_B7,"denv2"), "source": "PDB 4UT6 contact analysis <5A (DENV-2 numbering)", "role": "feature"},
    "FLE_fusion_loop": {"residues": mp(FLE,"denv2"), "source": "canonical fusion loop 98-110", "role": "feature"},
    "glycan_sites": {"residues": mp(GLY,"denv2"), "source": "N67/N153 E glycosylation (EDE dependence)", "role": "feature"},
    "escape_2D22": {"residues": [323], "source": "de Alwis 2012 PNAS 1200566109 (R323G; DENV-2 native)", "role": "ground_truth"},
    "2D22_associated": {"residues": [282,362], "source": "de Alwis 2012 PNAS Fig 3H (H282, D362 within 2D22/CR4354 footprint region)", "role": "ground_truth"}
  }
}
json.dump(epitopes, open("results/epitopes.json","w"), indent=1)
for s in epitopes:
    for k,v in epitopes[s].items(): print(s, k, v["residues"][:20], v["role"])
