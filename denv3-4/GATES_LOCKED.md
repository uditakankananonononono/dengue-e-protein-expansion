# GATES_LOCKED.md - DENV-3 / DENV-4 E-protein expansion slice (builder 15)
Locked: 2026-09-23 18:50 IST (before any outcome data is touched)

## Specific contribution (named before gates lock)
An epitope-resolved mutation atlas for DENV-3 and DENV-4 envelope (E) protein across
global isolates collected 2023-present, plus a quantified scoring tool:
**EERS - Epitope Exposure Risk Score** - a per-substitution score combining
(1) observed frequency in recent isolates, (2) membership in curated antibody epitopes
(from published mAb escape maps + PDB antibody-E structures), (3) cross-serotype
conservation context. Benchmark: EERS top-N enrichment of published escape mutations
vs two baselines (frequency-only ranking; conservation-only ranking), reported as
enrichment factor and AUROC. Prior art benchmarked against: Katzelnick et al.
antigenic maps (PMID-linked), Nextstrain dengue builds, published mAb escape studies
(5J7/DENV-3; 5H2/DENV-4; EDE/FLE cross-serotype).

## Success gates
G1 (data): >=200 DENV-3 and >=100 DENV-4 near-full-length or E-gene sequences with
   collection date >=2023-01-01 and country metadata, from NCBI Virus. If counts fall
   short, report honest counts and narrow claims; no padding.
G2 (positive control): pipeline recovers >=3 published escape/epitope mutations per
   serotype from the curated epitope set (e.g., DENV-4 mAb 5H2 escape K174; DENV-3
   5J7 epitope cluster; EDE N153 glycan dependence; FLE W101) - i.e., these residues
   rank high under EERS for mechanistic reasons, and known variable residues in the
   literature appear as genuinely variable in our alignment.
G3 (atlas): per-residue variability map of E (1-495) for both serotypes, with 2023+
   substitution frequencies, country spread, and epitope annotation. Real figures.
G4 (benchmark): EERS vs baselines reported with AUROC + top-20 enrichment on published
   escape sites. If EERS does not beat baselines, that negative result is reported as-is.
G5 (cross-serotype synthesis): comparison of expansion patterns DENV-3 vs DENV-4
   (which domains/epitopes accumulate substitutions), tied to reinfection/ADE
   mechanistic hypotheses with literature citations.

## Hard rules
- Real public data only (NCBI Virus, Nextstrain, RCSB PDB); source URLs + checksums in
  data manifest. No stubs, no simulations as results, no placeholder figures.
- Negative results preserved and reported; no re-fishing after seeing outcomes.
- No money spent; nothing sent as the user; compute stays in sandbox.
