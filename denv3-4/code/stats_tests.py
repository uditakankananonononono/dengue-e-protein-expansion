#!/usr/bin/env python3
import json, csv
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr
WORK="/home/sandbox/work/denv3-4/results"
def load(s):
    rows=list(csv.DictReader(open(f"{WORK}/atlas_{s}_v2.csv")))
    return [{**r,"pos":int(r["pos"]),"freq":float(r["freq"]),"epitopes":int(r["epitopes"])} for r in rows]
out={}
prof={}
for s in ("denv3","denv4"):
    t=load(s)
    diii=[r["freq"] for r in t if r["domain"]=="DIII"]
    rest=[r["freq"] for r in t if r["domain"]!="DIII"]
    u,p=mannwhitneyu(diii,rest,alternative="greater")
    # enrichment of ground-truth recovery among epitope-annotated vs not
    gt=json.load(open(f"{WORK}/epitopes.json"))[s]
    gtp={pp for k,v in gt.items() if v["role"]=="ground_truth" for pp in v["residues"]}
    variable={r["pos"] for r in t if r["freq"]>0.005}
    a=len(gtp & variable); b=len(gtp - variable)
    epipos={pp for k,v in gt.items() if v["role"]=="feature" for pp in v["residues"]}
    c=len((epipos-gtp) & variable); d=len((epipos-gtp) - variable)
    odds,pf=fisher_exact([[a,b],[c,d]])
    out[s]={"DIII_vs_rest_MannWhitney_p":p,"DIII_mean":sum(diii)/len(diii),"rest_mean":sum(rest)/len(rest),
            "GT_in_variable":a,"GT_total":a+b,"epitope_recovery_fisher_p":pf,"epitope_recovery_odds":odds if odds!=float('inf') else "inf"}
    prof[s]={r["pos"]:r["freq"] for r in t}
# cross-serotype profile correlation on mapped shared positions
maps=json.load(open("/tmp/e_posmaps.json"))
xs=[]; ys=[]
inv3={v:int(k) for k,v in maps["denv3"].items()}; inv4={v:int(k) for k,v in maps["denv4"].items()}
shared=set(inv3)&set(inv4)
for p in shared:
    xs.append(prof["denv3"].get(p,0)); ys.append(prof["denv4"].get(p,0))
rho,pr=spearmanr(xs,ys)
out["cross_serotype_profile_spearman"]={"rho":rho,"p":pr,"n":len(xs)}
json.dump(out, open(f"{WORK}/stats.json","w"), indent=1)
print(json.dumps(out,indent=1))
