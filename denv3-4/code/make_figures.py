#!/usr/bin/env python3
import json, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
WORK="/home/sandbox/work/denv3-4/results"
figd=f"{WORK}/figures"
summ=json.load(open(f"{WORK}/atlas_v2_summary.json"))
bench=json.load(open(f"{WORK}/benchmark.json"))
cohort=json.load(open(f"{WORK}/cohort_summary.json"))
ep=json.load(open(f"{WORK}/epitopes.json"))
DOM_BOUNDS=[("DI",1,52),(None,53,132),("DII",53,132),(None,133,193),("DI",133,193),(None,194,280),("DII",194,280),(None,281,297),("DI",281,297),("DIII",298,394),("stem",395,495)]
DOMS=[("DI",[(1,52),(133,193),(281,297)],"#dbe9f6"),("DII",[(53,132),(194,280)],"#fde9d9"),("DIII",[(298,394)],"#e2f0d9"),("stem/TM",[(395,495)],"#eeeeee")]
def load(s):
    rows=list(csv.DictReader(open(f"{WORK}/atlas_{s}_v2.csv")))
    return [{**r,"pos":int(r["pos"]),"freq":float(r["freq"]),"EERS":float(r["EERS"]),
             "epitopes":int(r["epitopes"]),"drift_vs_reference":r["drift_vs_reference"]=="True"} for r in rows]
gt=json.load(open(f"{WORK}/atlas_v2_summary.json"))
for s,label in [("denv3","DENV-3 (N=955, 2023+)"),("denv4","DENV-4 (N=129, 2023+)")]:
    t=load(s); L=len(t)
    fig,ax=plt.subplots(figsize=(11,4.2))
    for d,ranges,col in DOMS:
        for a,b in ranges: ax.axvspan(a,b,color=col,zorder=0)
    xs=[r["pos"] for r in t]; ys=[r["freq"] for r in t]
    ax.bar(xs,ys,width=1.0,color="#2c5d8a",zorder=3)
    for r in t:
        if r["drift_vs_reference"]:
            ax.plot(r["pos"],r["freq"],marker="v",color="purple",markersize=9,zorder=5)
    gts=[r for r in t if r["pos"] in {p for k,v in ep[s].items() if v["role"]=="ground_truth" for p in v["residues"]} and r["freq"]>0.005]
    for r in gts:
        ax.plot(r["pos"],r["freq"],marker="*",color="red",markersize=13,zorder=6)
        ax.annotate(str(r["pos"]),(r["pos"],r["freq"]),textcoords="offset points",xytext=(0,7),ha="center",fontsize=8,color="red")
    for d,ranges,col in DOMS:
        for a,b in ranges: ax.text((a+b)/2,-0.13*max(ys),d,ha="center",fontsize=9,color="#555")
    ax.set_xlim(0,L+1); ax.set_ylim(bottom=0)
    ax.set_xlabel("E protein position"); ax.set_ylabel("fraction of isolates with minority allele")
    ax.set_title(f"Per-residue variability of E protein, {label}\n(red star = published antibody-escape/natural-variation site; purple triangle = lineage-replacement drift vs reference)")
    plt.tight_layout(); plt.savefig(f"{figd}/F_variability_{s}.png",dpi=150); plt.close()
# F3 cohort
fig,axes=plt.subplots(2,2,figsize=(11,7))
for i,(s,ttl) in enumerate([("denv3","DENV-3"),("denv4","DENV-4")]):
    g=cohort[s]["geos"]; c=cohort[s]["clades"]
    ax=axes[i][0]; ks=list(g.keys())[:12]; ax.barh(ks[::-1],[g[k] for k in ks][::-1],color="#2c5d8a")
    ax.set_title(f"{ttl}: isolates by geography"); ax.tick_params(labelsize=7)
    ax=axes[i][1]; ks=list(c.keys())[:10]; ax.barh(ks[::-1],[c[k] for k in ks][::-1],color="#7a5195")
    ax.set_title(f"{ttl}: isolates by clade (nextclade)"); ax.tick_params(labelsize=7)
plt.tight_layout(); plt.savefig(f"{figd}/F_cohort.png",dpi=150); plt.close()
# F4 benchmark (honest)
fig,ax=plt.subplots(figsize=(7,4))
labels=["EERS","frequency-only","conservation-only"]
x=np.arange(3); w=0.35
d3=[bench["denv3"]["AUROC_EERS"],bench["denv3"]["AUROC_freq_only"],bench["denv3"]["AUROC_cons_only"]]
d4=[bench["denv4"]["AUROC_EERS"],bench["denv4"]["AUROC_freq_only"],bench["denv4"]["AUROC_cons_only"]]
ax.bar(x-w/2,d3,w,label="DENV-3",color="#2c5d8a"); ax.bar(x+w/2,d4,w,label="DENV-4",color="#c55a11")
ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylabel("AUROC (recover published escape/variable sites)")
ax.axhline(0.5,ls="--",color="grey")
for xi,v in zip(x,d3): ax.text(xi-w/2,v+0.01,f"{v:.2f}",ha="center",fontsize=8)
for xi,v in zip(x,d4): ax.text(xi+w/2,v+0.01,f"{v:.2f}",ha="center",fontsize=8)
ax.set_title("Benchmark: EERS vs baselines (honest result - no AUROC gain over frequency)")
ax.legend(); plt.tight_layout(); plt.savefig(f"{figd}/F_benchmark.png",dpi=150); plt.close()
# F5 epitope-group mean variability cross-serotype
fig,ax=plt.subplots(figsize=(8,4))
sets=["5J7_quaternary","EDE1_C8","EDE2_B7","FLE_fusion_loop","glycan_sites"]
d3=[summ["denv3"]["epitope_mean_freq"].get(k,0) for k in sets]
d4=[summ["denv4"]["epitope_mean_freq"].get(k,0) for k in sets]
x=np.arange(len(sets))
ax.bar(x-0.2,d3,0.4,label="DENV-3",color="#2c5d8a"); ax.bar(x+0.2,d4,0.4,label="DENV-4",color="#c55a11")
ax.set_xticks(x); ax.set_xticklabels(["5J7\n(DENV-3 specific)","EDE1","EDE2","fusion loop\n(FLE)","glycan\nsites"],fontsize=8)
ax.set_ylabel("mean minority-allele fraction per residue")
ax.set_title("Broadly cross-reactive epitopes (EDE/FLE/glycan) are frozen while\ntype-specific and DIII regions accumulate substitutions")
ax.legend(); plt.tight_layout(); plt.savefig(f"{figd}/F_epitope_groups.png",dpi=150); plt.close()
# F6 headline alleles
fig,axes=plt.subplots(1,2,figsize=(9,4))
gv3=gt["denv3"]["gt_verdicts"]; gv4=gt["denv4"]["gt_verdicts"]
ax=axes[0]
v=gv3["wahala_natural_variation:380"]; alleles={v["majority"]:v["maj_freq"]}; alleles.update({k:n/955 for k,n in v["subs"].items()})
ax.bar(alleles.keys(),alleles.values(),color="#2c5d8a")
ax.set_title("DENV-3 E380 (8A1 lateral-ridge site)\nI allele at 17.6% of 2023+ isolates"); ax.set_ylabel("allele fraction")
ax=axes[1]
v=gv4["escape_5H2:174"]; alleles={v["majority"]:v["maj_freq"]}; alleles.update({k:n/129 for k,n in v["subs"].items()})
ax.bar(alleles.keys(),alleles.values(),color="#c55a11")
ax.set_title("DENV-4 E174 (5H2 escape site)\nescape allele E is now the plurality state")
plt.tight_layout(); plt.savefig(f"{figd}/F_headline_alleles.png",dpi=150); plt.close()
print("figures written")
import os
print(os.listdir(figd))
