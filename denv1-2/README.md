# denv1-2: DENV-1/DENV-2 E-protein mutation atlas + EERS (builder 14 slice)

Reproduce end-to-end (<15 min, 2-core/2GB):
1. Data: NCBI Datasets CLI 18.37.0, taxids 11053/11060, collection_date>=2023-01-01,
   length>=1400nt -> data/ (byte-locked, see data/DATA_MANIFEST.sha256).
2. Nextclade 3.23.0, datasets community/v-gen-lab/dengue/denv1 + denv2:
   alignment, clades, aa substitution calls.
3. python3 code/pdb_contacts.py (structural epitope footprints from PDB 4C2I/4UIF) &&
   code/fix_mapping.py && code/build_epitopes.py && code/atlas_analysis.py &&
   code/atlas_v2.py && code/stats_tests.py && code/make_figures.py && code/make_tables.py
   (python 3.10, biopython/pandas/numpy/scipy/matplotlib)
4. Paper: pandoc paper/paper.md --pdf-engine=pdflatex -> paper/paper.pdf.
Gates: GATES_LOCKED.md (standalone first commit, before any data was fetched).
Results manifest: results/RESULTS_MANIFEST.sha256.
Findings: DENV-2 E71A (2D22/EDE-footprint residue) is fixed in the dominant 2II_F.1.1 lineage
(91.8% of all 1577 QC-passed 2023+ isolates; 99.4% majority frequency among the 872 isolates
with complete aligned E translations);
DENV-1 E155 (1F4 footprint) minority allele at 37%; DENV-1 E329 T329A escape allele
(PMC3312934) observed at low frequency in the wild. Cross-reactive EDE/FLE epitopes frozen.
Honest negative: EERS shows no AUROC gain over baselines; published in-vitro escape sites
are largely conserved in 2023+ field isolates (escape is rare in nature, not ranking failure).
