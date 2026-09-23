#!/usr/bin/env python3
"""DENV-3/4 E-protein mutation atlas + EERS benchmark. Reads nextclade outputs + NCBI metadata."""
import json, csv, re
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Align import PairwiseAligner

WORK="/home/sandbox/work/denv3-4"
# ---- per-serotype E consensus (from nextclade aligned translations) ----
def consensus(fasta):
    recs=list(SeqIO.parse(fasta,"fasta")); L=len(recs[0].seq); out=[]
    for i in range(L):
        c=Counter(str(r.seq[i]) for r in recs)
        for bad in ('-','X'): c.pop(bad,None)
        out.append(c.most_common(1)[0][0] if c else 'X')
    return "".join(out)
cons={"denv3":consensus("/tmp/ncout3/nextclade.cds_translation.E.fasta"),
      "denv4":consensus("/tmp/ncout4/nextclade.cds_translation.E.fasta")}
maps=json.load(open("/tmp/e_posmaps.json"))  # d2pos(str)->serotype pos
inv={s:{v:int(k) for k,v in m.items()} for s,m in maps.items()}  # serotype pos -> d2 pos

# ---- epitope sets (features use STRUCTURAL sets only; escape sets = ground truth) ----
EDE2_B7=[68,69,70,71,72,73,74,82,97,98,99,101,102,103,104,105,113,152,153,154,155,156,157,245,246,247,248,249]
EDE1_C8=[68,69,70,71,72,73,74,77,83,84,97,98,99,100,101,102,103,104,105,106,113,115,148,158,246,247,248,249,274,275,309,310,311,323,362]
FLE=list(range(98,111)); GLY=[67,153]
J7=[50,51,52,53,54,55,58,73,74,101,106,123,126,128,130,131,133,134,148,196,198,200,201,223,224,227,274,276,307,308,309]
def mp(l,s): return sorted({maps[s][str(p)] for p in l if str(p) in maps[s]})
SRC={"5J7_quaternary":"Fibriansah 2015 Science PMC4346626 Table 1",
     "EDE1_C8":"PDB 4UTA contact analysis <5A (this work)",
     "EDE2_B7":"PDB 4UT6 contact analysis <5A (this work)",
     "FLE_fusion_loop":"canonical fusion loop 98-110 (DENV-2 numbering mapped)",
     "glycan_sites":"N67/N153 glycosylation sites (EDE dependence)"}
EP={}
for s in ("denv3","denv4"):
    EP[s]={"EDE1_C8":mp(EDE1_C8,s),"EDE2_B7":mp(EDE2_B7,s),"FLE_fusion_loop":mp(FLE,s),"glycan_sites":mp(GLY,s)}
EP["denv3"]["5J7_quaternary"]=J7  # native DENV-3 numbering
GT={"denv3":{"escape_5J7_insertion":[269,270],
             "wahala_natural_variation":[301,302,329,380,386],
             "FL_variant_2023":[103,106]},
    "denv4":{"escape_5H2":[174],"escape_DV4-E75":[330,361,364]}}
ep_out={s:{**{k:{"residues":v,"source":SRC[k],"role":"feature"} for k,v in EP[s].items()},
           **{k:{"residues":v,"source":gtsrc,"role":"ground_truth"} for k,(gtsrc,v) in {
               "denv3":{"escape_5J7_insertion":("de Alwis 2012 PNAS 1200566109",[269,270]),
                        "wahala_natural_variation":("Wahala 2010 PLoS Pathog PMC2841629",[301,302,329,380,386]),
                        "FL_variant_2023":("eLife 87555 PMC10508882",[103,106])},
               "denv4":{"escape_5H2":("JVI 2007 5H2 escape K174",[174]),
                        "escape_DV4-E75":("Sukupolvi-Petty 2013 JVI PMC3754038",[330,361,364])}}[s].items()}} for s in ("denv3","denv4")}
json.dump(ep_out, open(f"{WORK}/results/epitopes.json","w"), indent=1)

# ---- cross-serotype conservation: align consensuses ----
al=PairwiseAligner(); al.mode="global"; al.match_score=2; al.mismatch_score=-1
al.open_gap_score=-8; al.extend_gap_score=-0.5
aln=al.align(cons["denv3"],cons["denv4"])[0]
xcons={}
for (a0,a1),(b0,b1) in zip(*aln.aligned):
    for i in range(int(a1)-int(a0)):
        xcons[("denv3",int(a0)+i+1)] = cons["denv3"][int(a0)+i]==cons["denv4"][int(b0)+i]
        xcons[("denv4",int(b0)+i+1)] = cons["denv3"][int(a0)+i]==cons["denv4"][int(b0)+i]

# ---- metadata ----
meta={}
for s,f in [("denv3",f"{WORK}/data/denv3_ncbi_report.jsonl"),("denv4",f"{WORK}/data/denv4_ncbi_report.jsonl")]:
    for line in open(f):
        r=json.loads(line)
        meta[r["accession"].split(".")[0]]={"geo":r.get("location",{}).get("geographic_location","?"),
          "region":r.get("location",{}).get("geographic_region","?"),
          "date":r.get("isolate",{}).get("collection_date",""), "len":r.get("length",0)}

# ---- per-position variability from nextclade aaSubstitutions ----
results={}
for s,tsv,n in [("denv3","/tmp/ncout3/nextclade.tsv",493),("denv4","/tmp/ncout4/nextclade.tsv",495)]:
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
    struct_sets=EP[s]
    epi_count={p:sum(1 for k,v in struct_sets.items() if p in v) for p in range(1,L+1)}
    gt_pos={p for v in GT[s].values() for p in v}
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
            "ground_truth":p in gt_pos})
    results[s]={"N":N,"L":L,"table":table,"clades":dict(clade_ct.most_common()),
                "geos":dict(geo_ct.most_common(25))}
    with open(f"{WORK}/results/atlas_{s}.csv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(table[0].keys())); w.writeheader()
        for row in table: w.writerow(row)

# ---- benchmark: AUROC + top-20 enrichment on ground-truth sites ----
def auroc(scores, labels):
    pairs=sorted(zip(scores,labels),key=lambda x:-x[0])
    P=sum(labels); Nn=len(labels)-P
    if P==0 or Nn==0: return None
    tp=fp=0; prev=None; area=0.0; last_fpr=last_tpr=0.0
    for sc,lb in pairs:
        if lb: tp+=1
        else: fp+=1
        tpr=tp/P; fpr=fp/Nn
        area+=(fpr-last_fpr)*(tpr+last_tpr)/2
        last_fpr,last_tpr=fpr,tpr
    return round(area,4)
bench={}
for s in ("denv3","denv4"):
    t=[r for r in results[s]["table"] if r["freq"]>0]
    labels=[1 if r["ground_truth"] else 0 for r in t]
    eers=[r["EERS"] for r in t]; fq=[r["freq"] for r in t]
    cn=[1.5 if r["cross_serotype_conserved"] else 1.0 for r in t]
    top20=sorted(t,key=lambda r:-r["EERS"])[:20]
    gt_total=sum(labels)
    bench[s]={"n_variable_positions":len(t),"n_ground_truth":gt_total,
        "AUROC_EERS":auroc(eers,labels),"AUROC_freq_only":auroc(fq,labels),
        "AUROC_cons_only":auroc(cn,labels),
        "top20_EERS_ground_truth_hits":sum(1 for r in top20 if r["ground_truth"]),
        "top20_freq_ground_truth_hits":sum(1 for r in sorted(t,key=lambda r:-r["freq"])[:20] if r["ground_truth"])}
    # positive control checks
    def ref_at(p): return cons[s][p-1]
    pc={}
    for p in [67,101,153]:
        row=results[s]["table"][p-1]
        pc[f"pos{p}_{ref_at(p)}"]={"ref":ref_at(p),"freq":row["freq"],"n_var":row["n_var"]}
    bench[s]["positive_controls_conserved"]=pc
    gt_rows={p:results[s]["table"][p-1] for p in gt_pos_set} if False else None
json.dump(bench, open(f"{WORK}/results/benchmark.json","w"), indent=1)
json.dump({s:{"N":results[s]["N"],"L":results[s]["L"],"clades":results[s]["clades"],"geos":results[s]["geos"]} for s in results},
          open(f"{WORK}/results/cohort_summary.json","w"), indent=1)
print(json.dumps(bench, indent=1))
for s in ("denv3","denv4"):
    print(s,"N=",results[s]["N"],"top EERS rows:")
    for r in sorted(results[s]["table"],key=lambda r:-r["EERS"])[:12]:
        print("  ",r["pos"],r["ref"],"f=",r["freq"],"epi=",r["n_struct_epitopes"],"xcons=",r["cross_serotype_conserved"],"EERS=",r["EERS"],"GT=",r["ground_truth"],"subs=",list(r["subs"].items())[:3])
