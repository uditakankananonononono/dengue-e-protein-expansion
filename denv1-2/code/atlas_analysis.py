#!/usr/bin/env python3
"""DENV-1/2 E-protein mutation atlas + EERS benchmark. Reads nextclade outputs + NCBI metadata."""
import json, csv, re
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Align import PairwiseAligner

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
maps=json.load(open("/tmp/e_posmaps_d12.json"))

ep_in=json.load(open(f"{WORK}/results/epitopes.json"))
EP={s:{k:v["residues"] for k,v in ep_in[s].items() if v["role"]=="feature"} for s in ep_in}
GT={s:{p for k,v in ep_in[s].items() if v["role"]=="ground_truth" for p in v["residues"]} for s in ep_in}

al=PairwiseAligner(); al.mode="global"; al.match_score=2; al.mismatch_score=-1
al.open_gap_score=-8; al.extend_gap_score=-0.5
aln=al.align(cons["denv1"],cons["denv2"])[0]
xcons={}
for (a0,a1),(b0,b1) in zip(*aln.aligned):
    for i in range(int(a1)-int(a0)):
        xcons[("denv1",int(a0)+i+1)] = cons["denv1"][int(a0)+i]==cons["denv2"][int(b0)+i]
        xcons[("denv2",int(b0)+i+1)] = cons["denv1"][int(a0)+i]==cons["denv2"][int(b0)+i]

meta={}
for s,f in [("denv1",f"{WORK}/data/denv1_ncbi_report.jsonl"),("denv2",f"{WORK}/data/denv2_ncbi_report.jsonl")]:
    for line in open(f):
        r=json.loads(line)
        meta[r["accession"].split(".")[0]]={"geo":r.get("location",{}).get("geographic_location","?"),
          "region":r.get("location",{}).get("geographic_region","?"),
          "date":r.get("isolate",{}).get("collection_date",""), "len":r.get("length",0)}

results={}
for s,tsv in [("denv1","/tmp/ncout1/nextclade.tsv"),("denv2","/tmp/ncout2/nextclade.tsv")]:
    rows=list(csv.DictReader(open(tsv),delimiter="\t"))
    rows=[r for r in rows if r.get("qc.overallStatus") in ("good","mediocre")]
    pos_subs=defaultdict(Counter); pos_geo=defaultdict(set); pos_clade=defaultdict(set)
    clade_ct=Counter(); geo_ct=Counter()
    for r in rows:
        acc=r["seqName"].split(".")[0]; m=meta.get(acc,{})
        clade_ct[r.get("clade","?")]+=1; geo_ct[m.get("geo","?")]+=1
        for sub in (r.get("aaSubstitutions") or "").split(","):
            mm=re.match(r'E:([A-Z])(\d+)([A-Z])$',sub.strip())
            if mm:
                p=int(mm.group(2))
                pos_subs[p][mm.group(3)]+=1
                pos_geo[p].add(m.get("geo","?")); pos_clade[p].add(r.get("clade","?"))
    N=len(rows); L=len(cons[s])
    epi_count={p:sum(1 for k,v in EP[s].items() if p in v) for p in range(1,L+1)}
    table=[]
    for p in range(1,L+1):
        nvar=sum(pos_subs[p].values()); f=nvar/N
        cons_flag=xcons.get((s,p),True)
        epi=epi_count[p]
        eers=f*(1+0.5*epi)*(1.5 if cons_flag else 1.0)
        table.append({"pos":p,"ref":cons[s][p-1],"n_var":nvar,"freq":round(f,5),
            "subs":dict(pos_subs[p].most_common(8)),"n_geo":len(pos_geo[p]),
            "n_clades":len(pos_clade[p]),"n_struct_epitopes":epi,
            "cross_serotype_conserved":cons_flag,"EERS":round(eers,6),
            "ground_truth":p in GT[s]})
    results[s]={"N":N,"L":L,"table":table,"clades":dict(clade_ct.most_common()),
                "geos":dict(geo_ct.most_common(25))}
    with open(f"{WORK}/results/atlas_{s}.csv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(table[0].keys())); w.writeheader()
        for row in table: w.writerow(row)

def auroc(scores, labels):
    pairs=sorted(zip(scores,labels),key=lambda x:-x[0])
    P=sum(labels); Nn=len(labels)-P
    if P==0 or Nn==0: return None
    tp=fp=0; area=0.0; last_fpr=last_tpr=0.0
    for sc,lb in pairs:
        if lb: tp+=1
        else: fp+=1
        tpr=tp/P; fpr=fp/Nn
        area+=(fpr-last_fpr)*(tpr+last_tpr)/2
        last_fpr,last_tpr=fpr,tpr
    return round(area,4)
bench={}
for s in ("denv1","denv2"):
    t=[r for r in results[s]["table"] if r["freq"]>0]
    labels=[1 if r["ground_truth"] else 0 for r in t]
    eers=[r["EERS"] for r in t]; fq=[r["freq"] for r in t]
    cn=[1.5 if r["cross_serotype_conserved"] else 1.0 for r in t]
    top20=sorted(t,key=lambda r:-r["EERS"])[:20]
    bench[s]={"n_variable_positions":len(t),"n_ground_truth":sum(labels),
        "AUROC_EERS":auroc(eers,labels),"AUROC_freq_only":auroc(fq,labels),
        "AUROC_cons_only":auroc(cn,labels),
        "top20_EERS_ground_truth_hits":sum(1 for r in top20 if r["ground_truth"]),
        "top20_freq_ground_truth_hits":sum(1 for r in sorted(t,key=lambda r:-r["freq"])[:20] if r["ground_truth"])}
    pc={}
    for p in [67,101,153]:
        row=results[s]["table"][p-1]
        pc[f"pos{p}_{row['ref']}"]={"ref":row["ref"],"freq":row["freq"],"n_var":row["n_var"]}
    bench[s]["positive_controls_conserved"]=pc
json.dump(bench, open(f"{WORK}/results/benchmark.json","w"), indent=1)
json.dump({s:{"N":results[s]["N"],"L":results[s]["L"],"clades":results[s]["clades"],"geos":results[s]["geos"]} for s in results},
          open(f"{WORK}/results/cohort_summary.json","w"), indent=1)
print(json.dumps(bench, indent=1))
for s in ("denv1","denv2"):
    print(s,"N=",results[s]["N"],"top EERS rows:")
    for r in sorted(results[s]["table"],key=lambda r:-r["EERS"])[:12]:
        print("  ",r["pos"],r["ref"],"f=",r["freq"],"epi=",r["n_struct_epitopes"],"xcons=",r["cross_serotype_conserved"],"EERS=",r["EERS"],"GT=",r["ground_truth"],"subs=",list(r["subs"].items())[:3])
