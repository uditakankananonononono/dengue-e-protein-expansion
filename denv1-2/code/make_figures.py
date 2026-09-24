#!/usr/bin/env python3
import json, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
WORK="./results"
figd=f"{WORK}/figures"
summ=json.load(open(f"{WORK}/atlas_v2_summary.json"))
bench=json.load(open(f"{WORK}/benchmark.json"))
cohort=json.load(open(f"{WORK}/cohort_summary.json"))
ep=json.load(open(f"{WORK}/epitopes.json"))
DOMS=[("DI",[(1,52),(133,193),(281,297)],"#dbe9f6"),("DII",[(53,132),(194,280)],"#fde9d9"),("DIII",[(298,394)],"#e2f0d9"),("stem/TM",[(395,495)],"#eeeeee")]
def load(s):
    rows=list(csv.DictReader(open(f"{WORK}/atlas_{s}_v2.csv")))
    return [{**r,"pos":int(r["pos"]),"freq":float(r["freq"]),"EERS":float(r["EERS"]),
             "epitopes":int(r["epitopes"]),"drift_vs_reference":r["drift_vs_reference"]=="True"} for r in rows]
N={"denv1":summ["denv1"]["N"],"denv2":summ["denv2"]["N"]}
for s,label in [("denv1",f"DENV-1 (N={N['denv1']}, 2023+)"),("denv2",f"DENV-2 (N={N['denv2']}, 2023+)")]:
    t=load(s); L=len(t)
    fig,ax=plt.subplots(figsize=(11,4.2))
    for d,ranges,col in DOMS:
        for a,b in ranges: ax.axvspan(a,b,color=col,zorder=0)
    xs=[r["pos"] for r in t]; ys=[r["freq"] for r in t]
    ax.bar(xs,ys,width=1.0,color="#2c5d8a",zorder=3)
    for r in t:
        if r["drift_vs_reference"]:
            ax.plot(r["pos"],r["freq"],marker="v",color="purple",markersize=9,zorder=5)
    gts=[r for r in t if r["pos"] in {p for k,v in ep[s].items() if v["role"]=="ground_truth" for p in v["residues"]} and r["freq"]>0.001]
    for r in gts:
        ax.plot(r["pos"],r["freq"],marker="*",color="red",markersize=13,zorder=6)
        ax.annotate(str(r["pos"]),(r["pos"],r["freq"]),textcoords="offset points",xytext=(0,7),ha="center",fontsize=8,color="red")
    for d,ranges,col in DOMS:
        for a,b in ranges: ax.text((a+b)/2,-0.13*max(ys),d,ha="center",fontsize=9,color="#555")
    ax.set_xlim(0,L+1); ax.set_ylim(bottom=0)
    ax.set_xlabel("E protein position"); ax.set_ylabel("fraction of isolates with minority allele")
    ax.set_title(f"Per-residue variability of E protein, {label}\n(red star = published antibody-escape/variable site; purple triangle = lineage-replacement drift vs reference)")
    plt.tight_layout(); plt.savefig(f"{figd}/F_variability_{s}.png",dpi=150); plt.close()
fig,axes=plt.subplots(2,2,figsize=(11,7))
for i,(s,ttl) in enumerate([("denv1","DENV-1"),("denv2","DENV-2")]):
    g=cohort[s]["geos"]; c=cohort[s]["clades"]
    ax=axes[i][0]; ks=list(g.keys())[:12]; ax.barh(ks[::-1],[g[k] for k in ks][::-1],color="#2c5d8a")
    ax.set_title(f"{ttl}: isolates by geography"); ax.tick_params(labelsize=7)
    ax=axes[i][1]; ks=list(c.keys())[:10]; ax.barh(ks[::-1],[c[k] for k in ks][::-1],color="#7a5195")
    ax.set_title(f"{ttl}: isolates by clade (nextclade)"); ax.tick_params(labelsize=7)
plt.tight_layout(); plt.savefig(f"{figd}/F_cohort.png",dpi=150); plt.close()
fig,ax=plt.subplots(figsize=(7,4))
labels=["EERS","frequency-only","conservation-only"]
x=np.arange(3); w=0.35
d1=[bench["denv1"]["AUROC_EERS"],bench["denv1"]["AUROC_freq_only"],bench["denv1"]["AUROC_cons_only"]]
d2=[bench["denv2"]["AUROC_EERS"],bench["denv2"]["AUROC_freq_only"],bench["denv2"]["AUROC_cons_only"]]
ax.bar(x-w/2,d1,w,label="DENV-1",color="#2c5d8a"); ax.bar(x+w/2,d2,w,label="DENV-2",color="#c55a11")
ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylabel("AUROC (recover published escape/variable sites)")
ax.axhline(0.5,ls="--",color="grey")
for xi,v in zip(x,d1): ax.text(xi-w/2,v+0.01,f"{v:.2f}",ha="center",fontsize=8)
for xi,v in zip(x,d2): ax.text(xi+w/2,v+0.01,f"{v:.2f}",ha="center",fontsize=8)
ax.set_title("Benchmark: EERS vs baselines\n(no AUROC gain; GT sites largely conserved)",fontsize=10)
ax.legend(); plt.tight_layout(); plt.savefig(f"{figd}/F_benchmark.png",dpi=150); plt.close()
fig,ax=plt.subplots(figsize=(8,4))
sets=[k for k in summ["denv1"]["epitope_mean_freq"]]
d1=[summ["denv1"]["epitope_mean_freq"].get(k,0) for k in sets]
d2=[summ["denv2"]["epitope_mean_freq"].get(k,0) for k in sets]
x=np.arange(len(sets))
ax.bar(x-0.2,d1,0.4,label="DENV-1",color="#2c5d8a"); ax.bar(x+0.2,d2,0.4,label="DENV-2",color="#c55a11")
ax.set_xticks(x); ax.set_xticklabels([k.replace("_","\n") for k in sets],fontsize=7)
ax.set_ylabel("mean minority-allele fraction per residue")
ax.set_title("Epitope-group mean variability: type-specific vs cross-reactive sets")
ax.legend(); plt.tight_layout(); plt.savefig(f"{figd}/F_epitope_groups.png",dpi=150); plt.close()
fig,axes=plt.subplots(1,2,figsize=(9,4))
gv1=summ["denv1"]["gt_verdicts"]
ax=axes[0]
v=gv1["escape_E106:329"]; alleles={v["majority"]:v["maj_freq"]}; alleles.update({k:n/N["denv1"] for k,n in v["subs"].items()})
ax.bar(alleles.keys(),alleles.values(),color="#2c5d8a")
ax.set_title("DENV-1 E329 (E106 escape site):\nT329A escape allele present in 2023+ isolates",fontsize=9); ax.set_ylabel("allele fraction")
ax=axes[1]
gv2=summ["denv2"]["gt_verdicts"]
v71=summ["denv2"]["drift_sites"]
row71=[d for d in v71 if d[0]==71][0]
ax.bar([row71[1]+" (old ref)",row71[2]+" (2023+ majority)"],[1-row71[3],row71[3]],color="#c55a11")
ax.set_title(f"DENV-2 E71 (2D22/EDE footprint):\nE71A fixed in dominant 2II_F.1.1 lineage\n({row71[3]*100:.1f}% of 2023+ isolates)",fontsize=9)
plt.tight_layout(); plt.savefig(f"{figd}/F_headline_alleles.png",dpi=150); plt.close()
print("figures written")
import os
print(os.listdir(figd))
