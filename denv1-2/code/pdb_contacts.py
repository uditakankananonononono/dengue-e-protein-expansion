#!/usr/bin/env python3
"""All-atom <=6A contact analysis: Fab chains vs E-protein chains.
4UIF: E = A/C/E, Fab heavy = G/I/K, light = H/J/L (DENV-2 + 2D22).
4C2I: E = A/C/E, Fab heavy = H/M, light = L/N (DENV-1 + 1F4).
Outputs union of E residues in contact with any Fab chain."""
from collections import defaultdict
def atoms(f, chains):
    out = []
    for line in open(f):
        if line.startswith("ATOM") and line[21] in chains:
            try:
                out.append((int(line[22:26]), float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
    return out
for pdb, echains, fabchains in [("/tmp/4UIF.pdb", set("ACE"), set("GIKHJL")),
                                ("/tmp/4C2I.pdb", set("ACE"), set("HMLN"))]:
    E = atoms(pdb, echains); F = atoms(pdb, fabchains)
    cell = 7.0; grid = defaultdict(list)
    for (r,x,y,z) in E: grid[(int(x//cell),int(y//cell),int(z//cell))].append((r,x,y,z))
    contacts = set()
    for (r,x,y,z) in F:
        k = (int(x//cell),int(y//cell),int(z//cell))
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for (r2,x2,y2,z2) in grid.get((k[0]+dx,k[1]+dy,k[2]+dz),[]):
                        if (x-x2)**2+(y-y2)**2+(z-z2)**2 <= 36: contacts.add(r2)
    print(pdb, sorted(contacts))
