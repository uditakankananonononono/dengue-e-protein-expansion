# GATES_LOCKED.md - DENV-1 / DENV-2 E-protein expansion slice (builder 14)
Locked: 2026-09-24 09:37 IST (before any outcome data is touched)

## Specific contribution (named before gates lock)
An epitope-resolved mutation atlas for DENV-1 and DENV-2 envelope (E) protein across
global isolates collected 2023-present, plus the same quantified scoring tool used in the
sealed denv3-4 slice:
**EERS - Epitope Exposure Risk Score** - a per-substitution score combining
(1) observed frequency in recent isolates, (2) membership in curated antibody epitopes
(from published mAb escape maps + PDB antibody-E structures), (3) cross-serotype
conservation context. Benchmark: EERS top-N enrichment of published escape mutations
vs two baselines (frequency-only ranking; conservation-only ranking), reported as
enrichment factor and AUROC. Prior art benchmarked against: Katzelnick et al.
antigenic maps (PMID-linked), Nextstrain dengue builds, published mAb escape studies
(DENV-2 2D22 DIII lateral-ridge quaternary mAb; DENV-1 1F4 DI hinge mAb; DENV-1 DIII
neutralizing escape maps; EDE/FLE cross-serotype).

## Success gates
G1 (data): >=300 DENV-1 and >=400 DENV-2 near-full-length or E-gene sequences with
   collection date >=2023-01-01 and country metadata, from NCBI Virus. If counts fall
   short, report honest counts and narrow claims; no padding.
G2 (positive control): pipeline recovers >=3 published escape/epitope residues per
   serotype from the curated epitope set (e.g., DENV-2 2D22 epitope cluster;
   DENV-1 1F4 DI-hinge contact; EDE N153 glycan dependence; FLE W101) - i.e., these
   residues rank high under EERS for mechanistic reasons, and known variable residues
   in the literature appear as genuinely variable in our alignment.
G3 (atlas): per-residue variability map of E (1-495) for both serotypes, with 2023+
   substitution frequencies, country spread, and epitope annotation. Real figures.
G4 (benchmark): EERS vs baselines reported with AUROC + top-20 enrichment on published
   escape sites. If EERS does not beat baselines, that negative result is reported as-is.
G5 (cross-serotype synthesis): comparison of expansion patterns DENV-1 vs DENV-2
   (which domains/epitopes accumulate substitutions), tied to reinfection/ADE
   mechanistic hypotheses with literature citations.

## Hard rules
- Real public data only (NCBI Virus, Nextstrain, RCSB PDB); source URLs + checksums in
  data manifest. No stubs, no simulations as results, no placeholder figures.
- Negative results preserved and reported; no re-fishing after seeing outcomes.
- No money spent; nothing sent as the user; compute stays in sandbox.
