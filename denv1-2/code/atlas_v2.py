#!/usr/bin/env python3
"""Atlas v2: reference-drift vs within-population polymorphism, GT site verdicts,
substitution-level watchlist, domain/epitope expansion summary."""
import json, csv, re
from collections import Counter, defaultdict
from Bio import SeqIO

WORK="."
def consensus(fasta):
    recs=list(SeqIO.parse(fasta,"fasta")); L=max(len(r.seq) for r in recs); out=[]
    for i in range(L):
        c=Counter(str(r.seq[i]) for r in recs if len(r.seq)>i)
        c.pop('-',None); c.pop('X',None)
        out.append(c.most_common(1)[0][0] if c else 'X')
    return "".join(out)
cons={"denv1":consensus("/tmp/ncout1/nextclade.cds_translation.E.fasta"),
      "denv2":consensus("/tmp/ncout2/nextclade.cds_translation.E.fasta")}
ep=json.load(open(f"{WORK}/results/epitopes.json"))
feat={s:{k:v["residues"] for k,v in ep[s].items() if v["role"]=="feature"} for s in ep}
gt={s:{k:v["residues"] for k,v in ep[s].items() if v["role"]=="ground_truth"} for s in ep}
meta={}
for s,f in [("denv1",f"{WORK}/data/denv1_ncbi_report.jsonl"),("denv2",f"{WORK}/data/denv2_ncbi_report.jsonl")]:
    for line in open(f):
        r=json.loads(line); meta[r["accession"].split(".")[0]]=r.get("location",{}).get("geographic_location","?")
DOM=[("DI",[(1,52),(133,193),(281,297)]),("DII",[(53,132),(194,280)]),("DIII",[(298,394)]),("stem_anchor",[(395,495)])]
def domain(p):
    for d,rs in DOM:
        for a,b in rs:
            if a<=p<=b: return d
    return "?"
out={}
NC={"denv1":"/tmp/ncout1","denv2":"/tmp/ncout2"}
for s,tsv in [("denv1","/tmp/ncout1/nextclade.tsv"),("denv2","/tmp/ncout2/nextclade.tsv")]:
    rows=[r for r in csv.DictReader(open(tsv),delimiter="\t") if r.get("qc.overallStatus") in ("good","mediocre")]
    N=len(rows); L=len(cons[s])
    aligned={r.id.split('.')[0]:str(r.seq) for r in SeqIO.parse(f"{NC[s]}/nextclade.cds_translation.E.fasta","fasta")}
    keep=[r["seqName"].split(".")[0] for r in rows]
    pos_counts=[]
    for p in range(1,L+1):
        c=Counter(aligned[a][p-1] for a in keep if a in aligned and len(aligned[a])>=p and aligned[a][p-1] not in "-X")
        tot=sum(c.values())
        if not tot: c["X"]=0; tot=1
        maj,majn=c.most_common(1)[0]
        subs={aa:n for aa,n in c.items() if aa!=maj}
        pos_counts.append({"pos":p,"majority":maj,"maj_freq":majn/tot if tot else 0,
                           "subs":subs,"n_var":sum(subs.values()),"n_obs":tot})
    dsref={}
    for r in rows:
        for sub in (r.get("aaSubstitutions") or "").split(","):
            m=re.match(r'E:([A-Z])(\d+)([A-Z])$',sub.strip())
            if m: dsref.setdefault(int(m.group(2)),m.group(1))
    epi_count={p:sum(1 for v in feat[s].values() if p in v) for p in range(1,L+1)}
    table=[]
    for pc in pos_counts:
        p=pc["pos"]; f=pc["n_var"]/N
        drift = pc["maj_freq"]>=0.95 and p in dsref and dsref[p]!=pc["majority"]
        eers=f*(1+0.5*epi_count[p])
        table.append({**pc,"pos":p,"freq":round(f,5),"epitopes":epi_count[p],
                      "domain":domain(p),"drift_vs_reference":bool(drift),
                      "dataset_ref":dsref.get(p,pc["majority"]),"EERS":round(eers,6)})
    gtv={}
    for gname,plist in gt[s].items():
        for p in plist:
            pc=pos_counts[p-1]
            gtv[f"{gname}:{p}"]={"majority":pc["majority"],"maj_freq":round(pc["maj_freq"],4),
               "subs":pc["subs"],"dataset_ref":dsref.get(p,"?"),
               "verdict":("variable" if pc["n_var"]>=2 else ("drifted_vs_ref" if dsref.get(p,pc["majority"])!=pc["majority"] else "conserved"))}
    wl=[]
    for pc in pos_counts:
        p=pc["pos"]
        for aa,n in pc["subs"].items():
            wl.append({"pos":p,"from":pc["majority"],"to":aa,"n":n,"freq":round(n/N,5),
                       "epitopes":epi_count[p],"domain":domain(p),"EERS":round(n/N*(1+0.5*epi_count[p]),6)})
    wl.sort(key=lambda r:-r["EERS"])
    out[s]={"N":N,"gt_verdicts":gtv,"watchlist":wl[:40],
            "n_drift_sites":sum(1 for t in table if t["drift_vs_reference"]),
            "drift_sites":[(t["pos"],t["dataset_ref"],t["majority"],round(t["maj_freq"],3),t["epitopes"]) for t in table if t["drift_vs_reference"]]}
    with open(f"{WORK}/results/atlas_{s}_v2.csv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["pos","majority","maj_freq","n_var","n_obs","freq","epitopes","domain","drift_vs_reference","dataset_ref","EERS"])
        w.writeheader()
        for t in table:
            tt=dict(t); tt.pop("subs"); w.writerow(tt)
    dom_var=defaultdict(float); dom_n=defaultdict(int)
    for t in table:
        dom_var[t["domain"]]+=t["freq"]; dom_n[t["domain"]]+=1
    epi_var={k:round(sum(table[p-1]["freq"] for p in v)/len(v),5) for k,v in feat[s].items()}
    out[s]["domain_mean_freq"]={d:round(dom_var[d]/dom_n[d],5) for d in dom_n}
    out[s]["epitope_mean_freq"]=epi_var
json.dump(out, open(f"{WORK}/results/atlas_v2_summary.json","w"), indent=1)
for s in out:
    print("=====",s,"N=",out[s]["N"],"drift sites vs dataset-ref:",out[s]["n_drift_sites"])
    print("drift:",out[s]["drift_sites"])
    print("GT verdicts:")
    for k,v in out[s]["gt_verdicts"].items(): print("  ",k,v["verdict"],v["majority"],round(v["maj_freq"],3),"subs:",v["subs"],"dsref:",v["dataset_ref"])
    print("domain mean freq:",out[s]["domain_mean_freq"])
    print("epitope mean freq:",out[s]["epitope_mean_freq"])
    print("top watchlist:",[(w['pos'],w['from'],w['to'],w['n'],w['epitopes']) for w in out[s]['watchlist'][:10]])
