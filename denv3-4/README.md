# denv3-4: DENV-3/DENV-4 E-protein mutation atlas + EERS (builder 15 slice)

Reproduce end-to-end (<15 min, 2-core/2GB):
1. Data: NCBI Datasets CLI 18.37.0, taxids 11069/11070, collection_date>=2023-01-01,
   length>=1400nt -> data/ (byte-locked, see data/DATA_MANIFEST.sha256).
2. Nextclade 3.23.0, datasets community/v-gen-lab/dengue/denv3 + denv4 (2026-04-14 build):
   alignment, clades, aa substitution calls.
3. python3 code/build_epitopes.py && code/atlas_analysis.py && code/atlas_v2.py &&
   code/stats_tests.py && code/make_figures.py  (python 3.10, biopython/pandas/numpy/scipy/matplotlib)
4. Paper: pandoc paper/paper.md --pdf-engine=pdflatex -> paper/paper.pdf (19 pp).
Gates: GATES_LOCKED.md (standalone first commit). Results manifest: results/RESULTS_MANIFEST.sha256.
Findings: DENV-4 E174 5H2-escape allele plurality; DENV-3 E380 lateral-ridge polymorphism on
3III_B.3.2; cross-reactive epitopes frozen vs type-specific drift. Honest negative: EERS shows
no AUROC gain over frequency-only ranking (reported, not tuned).
