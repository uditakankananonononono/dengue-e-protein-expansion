#!/usr/bin/env python3
"""Contact-geometry analysis of DENV-2 E residue 71 (the E71A headline site) against the
antibody Fabs in PDB 4UIF (2D22), 4UTA (EDE1 C8) and 4UT6 (EDE2 B7).
For each structure: every Fab residue with an atom within 6A of any E71 atom, with the
minimum distance overall and the minimum distance to the Glu side-chain-only atoms
(CG, CD, OE1, OE2 - the atoms deleted by the E->A substitution). Classification:
 - close_sidechain_contact_lost_in_E71A: a deleted atom comes within 4A of the Fab
 - sidechain_contact_within_6A_weakened: deleted atoms within 6A but not <4A
 - backbone_or_CB_contact_retained: contact mediated only by atoms alanine keeps
Numbering anchors verified in all three structures: N67, W101, N153 (N153 unresolved in
4UTA chain A only). Outputs results/e71_structural.json."""
import json
from collections import defaultdict

STRUCTS = [
    ("/tmp/4UIF.pdb", set("ACE"), set("GIKHJL"), "2D22 (4UIF)"),
    ("/tmp/4UTA.pdb", set("AB"), set("HILM"), "EDE1 C8 (4UTA)"),
    ("/tmp/4UT6.pdb", set("AB"), set("HILM"), "EDE2 B7 (4UT6)"),
]
LOST = {"CG", "CD", "OE1", "OE2"}  # Glu atoms beyond CB, absent in Ala

def atoms(f, chains):
    out = []
    for line in open(f):
        if line.startswith("ATOM") and line[21] in chains:
            try:
                out.append(dict(chain=line[21], resnum=int(line[22:26]),
                                resname=line[17:20].strip(), atom=line[12:16].strip(),
                                x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54])))
            except ValueError:
                pass
    return out

report = {}
for pdb, echains, fabchains, label in STRUCTS:
    E = atoms(pdb, echains); F = atoms(pdb, fabchains)
    e71 = [a for a in E if a["resnum"] == 71]
    resname = e71[0]["resname"]
    cell = 7.0; grid = defaultdict(list)
    for a in F:
        grid[(int(a["x"]//cell), int(a["y"]//cell), int(a["z"]//cell))].append(a)
    hits = {}
    for a in e71:
        k = (int(a["x"]//cell), int(a["y"]//cell), int(a["z"]//cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for b in grid.get((k[0]+dx, k[1]+dy, k[2]+dz), []):
                        d = ((a["x"]-b["x"])**2+(a["y"]-b["y"])**2+(a["z"]-b["z"])**2)**0.5
                        if d <= 6.0:
                            key = (b["chain"], b["resnum"], b["resname"])
                            h = hits.setdefault(key, {"min_all": 9e9, "min_lost": 9e9, "lost_atoms": set()})
                            h["min_all"] = min(h["min_all"], d)
                            if a["atom"] in LOST:
                                h["min_lost"] = min(h["min_lost"], d); h["lost_atoms"].add(a["atom"])
    rows = []
    for (ch, rn, rname), h in sorted(hits.items(), key=lambda kv: kv[1]["min_all"]):
        cls = ("close_sidechain_contact_lost_in_E71A" if h["min_lost"] <= 4.0
               else "backbone_or_CB_contact_retained" if not h["lost_atoms"]
               else "sidechain_contact_within_6A_weakened")
        rows.append(dict(fab_chain=ch, fab_resnum=rn, fab_resname=rname,
                         min_dist_A=round(h["min_all"], 2),
                         min_dist_lost_atoms_A=(round(h["min_lost"], 2) if h["lost_atoms"] else None),
                         lost_e71_atoms=sorted(h["lost_atoms"]), classification=cls))
    report[label] = dict(pdb=pdb.split("/")[-1], e71_resname=resname,
                         fab_residues_within_6A=len(rows),
                         n_close_sidechain_lost=sum(1 for r in rows if r["classification"] == "close_sidechain_contact_lost_in_E71A"),
                         contacts=rows)

json.dump(report, open("results/e71_structural.json", "w"), indent=1)
for label, r in report.items():
    print(label, "E71:", r["e71_resname"], "| Fab res <6A:", r["fab_residues_within_6A"],
          "| close sidechain contacts lost in E71A:", r["n_close_sidechain_lost"])
