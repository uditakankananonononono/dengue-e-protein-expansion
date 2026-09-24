#!/usr/bin/env python3
"""Emit results markdown tables T1-T4 + TB_var from the analysis outputs."""
import json, csv
from collections import Counter
R="./results"
summ=json.load(open(f"{R}/atlas_v2_summary.json"))
cohort=json.load(open(f"{R}/cohort_summary.json"))
ep=json.load(open(f"{R}/epitopes.json"))
meta={}
for s,f in [("denv1","./data/denv1_ncbi_report.jsonl"),("denv2","./data/denv2_ncbi_report.jsonl")]:
    dates=[]; geos=set()
    for line in open(f):
        r=json.loads(line)
        dates.append(r.get("isolate",{}).get("collection_date",""))
        geos.add(r.get("location",{}).get("geographic_location","?").split(":")[0].strip())
    meta[s]={"dates":sorted(d for d in dates if d),"geos":geos}
acc={s:sum(1 for _ in open(f"./data/{s}_2023_accessions.txt")) for s in ("denv1","denv2")}
def domline(s):
    c=cohort[s]["clades"]
    k,v=max(c.items(), key=lambda kv: kv[1])
    return f"{k} ({v})"
with open(f"{R}/T1_cohort.md","w") as f:
    f.write("| Quantity | DENV-1 | DENV-2 |\n|---|---|---|\n")
    f.write(f"| Accessions downloaded (2023+, >=1400 nt) | {acc['denv1']} | {acc['denv2']} |\n")
    f.write(f"| QC-passed isolates in analysis | {summ['denv1']['N']} | {summ['denv2']['N']} |\n")
    f.write(f"| Collection-date span | 2023 to {meta['denv1']['dates'][-1]} | 2023 to {meta['denv2']['dates'][-1]} |\n")
    f.write(f"| Distinct countries | {len(meta['denv1']['geos'])} | {len(meta['denv2']['geos'])} |\n")
    f.write(f"| Dominant lineage | {domline('denv1')} | {domline('denv2')} |\n")
with open(f"{R}/T2_epitopes.md","w") as f:
    f.write("| Serotype | Epitope set | Role | # residues | Source |\n|---|---|---|---|---|\n")
    for s in ("denv1","denv2"):
        for k,v in ep[s].items():
            f.write(f"| {s.upper()} | {k} | {v['role']} | {len(v['residues'])} | {v['source']} |\n")
with open(f"{R}/T3_gt.md","w") as f:
    f.write("| Serotype | Ground-truth site (source) | Majority allele (freq) | Minority alleles | Verdict |\n|---|---|---|---|---|\n")
    src={s:{k:v["source"] for k,v in ep[s].items()} for s in ep}
    for s in ("denv1","denv2"):
        for gname,vv in summ[s]["gt_verdicts"].items():
            g,p=gname.rsplit(":",1)
            subs=", ".join(f"{k}:{n}" for k,n in vv["subs"].items()) or "-"
            f.write(f"| {s.upper()} | {g}:{p} ({src[s][g].split('(')[0].strip()}) | {vv['majority']} ({vv['maj_freq']}) | {subs} | {vv['verdict']} |\n")
for s in ("denv1","denv2"):
    with open(f"{R}/T4_watchlist_{s}.md","w") as f:
        f.write("| E pos | from | to | n isolates | fraction | #structural epitopes | domain | EERS |\n|---|---|---|---|---|---|---|---|\n")
        for w in summ[s]["watchlist"][:15]:
            f.write(f"| {w['pos']} | {w['from']} | {w['to']} | {w['n']} | {w['freq']} | {w['epitopes']} | {w['domain']} | {w['EERS']} |\n")
    rows=[{**r,"pos":int(r["pos"]),"freq":float(r["freq"]),"n_var":int(r["n_var"])} for r in csv.DictReader(open(f"{R}/atlas_{s}_v2.csv"))]
    rows=[r for r in rows if r["n_var"]>0]
    rows.sort(key=lambda r:-r["freq"])
    aligned=json.load(open(f"{R}/../results/epitopes.json"))  # placeholder
    # minority allele detail from watchlist
    wl={(w["pos"],w["to"]):w for w in summ[s]["watchlist"]}
    with open(f"{R}/TB_var_{s}.md","w") as f:
        f.write("| E pos | majority | minority alleles (n) | fraction | domain | #epitopes | drift vs ref |\n|---|---|---|---|---|---|---|\n")
        for r in rows[:15]:
            subs="; ".join(f"{w['to']}({w['n']})" for w in summ[s]["watchlist"] if w["pos"]==r["pos"])
            f.write(f"| {r['pos']} | {r['majority']} | {subs} | {r['freq']} | {r['domain']} | {r['epitopes']} | {r['drift_vs_reference']} |\n")
print("tables done")
