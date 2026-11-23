# Supplementary Data S1 — raw qPCR data, metadata and analysis code

This folder reproduces, from the original instrument files, every number, table and figure in the manuscript
*Mitochondrial DNA integrity and gene expression across the brooding period and among subclades of the Heliopora
species complex*.

```
Supplementary_Data_S1/
├── README.md                    this guide
├── requirements.txt             Python packages
├── run_all.py                   runs steps 1-9 in order
├── data/
│   ├── metadata/samples.csv     one row per colony: site, date, subclade, sex, sampling period, field data, DNA dilution, flags
│   └── raw/                     original files, renamed only (never edited)
│       ├── gene_expression/     LinRegPCR outputs of the two cDNA plates; plate map; eds/ (original run files)
│       ├── mtDNA/               SRP54, short COI and long COI plates (instrument exports, raw curves, LinRegPCR); eds/
│       ├── mtDNA_KR_replicate_20180314/   independent repeat run (KR colonies; ACTB, SRP54, short and long COI); eds/
│       └── sample_records/      primer list; mtDNA sample sheet (DNA concentrations and dilutions)
├── scripts/                     qpcr_lib.py (shared functions) + steps 01-09
└── results/                     created by run_all.py (tables .csv, figures .svg/.png, SUMMARY.md, all_results.xlsx)
```

## 1. Quick start

```
pip install -r requirements.txt
python run_all.py
```

Then open `results/SUMMARY.md`. The figures are `results/Fig1_expression.svg`, `results/Fig2_mtDNA.svg` and
`results/FigS1_copy_number_vs_damage.svg` (Supplementary Figure S1).
The run takes about one minute (tested with Python 3.13 and the package versions in `requirements.txt`).

## 2. Names used

- **Subclades** follow Fig. 1 of Taninaka et al. (2021): coe.3 = historical HC-B, sp.2 = HC-A′, sp.3 = HC-A.
  The sample sheets and the result tables use the historical labels B, A′ and A; `samples.csv` gives both
  (columns `lineage` and `subclade`).
- **Sampling periods.** "Before brooding" = 13 July 2014; "during brooding" = the population's reproductive period
  (21, 25 or 28 July 2014). Male colonies do not brood; the term refers to the period, not to individual colonies.
  In the result tables this period is labelled `reproductive`.
- **Assays.** Short COI (~100 bp) and long COI (1,124 bp) are mitochondrial; SRP54 and ACTB are nuclear references.

## 3. The study in one table

| Data set | Site | Before brooding (13 Jul 2014) | During brooding | Steps |
|---|---|---|---|---|
| Gene expression | Kohama Reef (KR) | SHC417 ♂, SHC418 ♀, SHC419 ♂ | SHC621 ♀, SHC622 ♂, SHC623 ♂ (25 Jul) | 1 |
| mtDNA | KR | SHC413 ♀, SHC414 ♀, SHC415 ♂ | SHC608, SHC612, SHC614 (21 Jul) | 2-4 |
| mtDNA | Yonara Iriomote Slope (YIS) | SHC422 ♂?, SHC423 ♀?, SHC424 | SHC593, SHC595, SHC596 (21 Jul) | 2-3 |
| mtDNA | Yonara Kohama Slope (YKS) | SHC426, SHC427, SHC428 | SHC699 ♀, SHC701 ♀, SHC706 ♀ (28 Jul) | 2-3 |
| mtDNA subclade panel | Yonara, 13 Jun 2017 | sp.3: SHC1170-72; sp.2: SHC1180-82; coe.3: SHC1160-62 | — | 2-3 |
| Exploratory (Table S2) | YKS, 17 Jul 2014 (sp.3, after gamete release); archived DNA 2007-2011 (sp.2, sp.3) | SHC518, SHC519, SHC524; AMM, KUM, OUR, ZIN, KAT | — | 3, 6 |
| Copy number v damage (Fig. S1) | the 27 colonies of the four mtDNA rows above (exploratory row not used) | — | — | 8 |
| coe.3, 2014 v 2017 | 2017 coe.3 (3) against 2014 coe.3 before brooding (9; and YIS before alone, 3) | — | — | 9 |

All 2014 colonies in the main comparisons are coe.3. Expression and mtDNA colonies are different individuals from the
same reef, subclade and week, so the two data sets are compared at the population level, not colony by colony.
`samples.csv` was compiled from the field sheets and sample lists of the sampling campaigns.

## 4. Principles used throughout

1. **The colony is the unit.** Technical duplicates (or triplicates) are averaged on the log2 scale before any test;
   counting wells as replicates would inflate n.
2. **Efficiency-corrected quantities.** `log2 N0 = log2(threshold) − Cq × log2(E)`, where E is the mean efficiency of
   the amplicon (LinRegPCR rule: mean of wells within 5 % of the median). The threshold term cancels in every ratio.
3. **Two-sided tests with the null hypothesis stated** in each script's header.
4. **Every decision is in code** (outlier rule, sample exclusions, relabelling), not made by hand.

## 5. Steps

### Step 1 — Gene expression (`scripts/01_expression.py`)

Does expression at KR differ between before (13 Jul) and during brooding (25 Jul)?

1. Reads the two LinRegPCR files. Plate 1: ACTB, RPL12, D2R, mtMutS, CRY1, COI. Plate 2: SRP54, CAT, GPX, HSP70,
   HSP90, SOD and ACTB again (same cDNA). Plate labels C1-C3 (before) and T1-T3 (during) map to colonies via `samples.csv`.
2. **Outlier rule.** One ACTB well (plate 1, C1, duplicate 1) differs by 1.1 cycles from its duplicate; the same cDNA
   on plate 2 agrees with the other duplicate, so this well is excluded. No other well meets the rule.
3. **Reference genes** (`01_reference_genes.csv`). ACTB + RPL12 are used. SRP54 is not used: efficiency 1.76, Cq about
   31, poor duplicate agreement, and drift relative to ACTB and RPL12.
4. **Normalisation.** For each colony: target log2 N0 minus the mean of ACTB (both plates) and RPL12.
5. **Test.** Student t (df = 4), 95 % CI, exact permutation p (smallest possible 0.10 with 3 v 3), and
   Benjamini-Hochberg q over the nine targets.
6. **Checks.** Every reference set, Welch t, individual-well efficiencies, outlier kept, E = 2, leave-one-colony-out,
   and colonies needed for 80 % power (`01_expression_*.csv`).

*Worked example, SOD.* From `01_expression_per_colony_log2NRQ.csv`: before = 1.39, 0.59, 0.08 (mean 0.69);
during = −0.71, −0.22, −1.08 (mean −0.67). Difference = −1.35 log2, fold change 2^−1.35 = 0.39.
Student t on these six numbers gives p = 0.041; with nine targets, q = 0.37.

### Step 2 — mtDNA quantities (`scripts/02_mtdna_quantities.py`)

Three assays on total DNA (equal input, 100 ng in 20 µl): SRP54 (nuclear), short COI and long COI.

1. **Validation.** LinRegPCR output exists for the long COI plate and for a first short-COI run (later replaced
   because of poor duplicate agreement). The curve-analysis code in `qpcr_lib.analyse_curves` is run on these two
   plates and compared with LinRegPCR (`02_estimator_validation.csv`): mean efficiency within 0.012, Cq identical up to
   a constant offset (SD 0.02-0.03).
2. **Quantities.** SRP54 (23 Mar 2018) and short COI (26 Mar 2018) from raw curves; long COI from LinRegPCR.
   Wells labelled "SHC520" are SHC524 (SHC520 is in no sample list).
3. **Diagnostics.** Duplicate agreement per assay (`02_mtdna_duplicate_agreement.csv`) and the sensitivity of SRP54
   Cq to the threshold setting (`02_srp54_threshold_sensitivity.csv`).

### Step 3 — Copy number and damage (`scripts/03_mtdna_statistics.py`)

- **Copy number:** `log2 copy number = log2 N0(short COI) − log2 N0(SRP54)`. SRP54-independent check: short COI per
  unit DNA input.
- **Damage** (semi-long run qPCR; Rothfuss et al. 2010): `LS = log2 N0(long COI) − log2 N0(short COI)`,
  `r = 2^(LS − mean LS of the reference group)`, lesions per 10 kb = `−ln(r) × 10,000 / 1,124` (Poisson) or
  `(1 − r) × 10,000 / 1,124` (Rothfuss). Reference group = before-brooding colonies of the same site, or coe.3 (2017).
- **Exclusions from copy number only** (SRP54 failed or amplified late): SHC701, AMM6, KAT9, KAT10, KAT12. They remain
  in the damage analysis, which does not use SRP54.
- **Tests:** Welch t within site; two-way ANOVA (site × period); KR vs YIS difference in change; one-way ANOVA and
  Tukey HSD for subclades.
- **Robustness:** leave-one-colony-out, each long-COI duplicate alone, E = 2.

### Step 4 — Independent KR repeat run (`scripts/04_kr_replicate_20180314.py`)

The 14 Mar 2018 run measured the six KR mtDNA colonies in triplicate with ACTB and SRP54. Copy number did not change
(ACTB-normalised 1.01-fold, SRP54-normalised 1.08-fold), and the two nuclear references agreed. Short-COI curves for
SHC614 are missing from the export, so damage in this run compares three with two colonies. No-template controls are
listed in `04_ntc.csv`.

### Step 5 — Summary (`scripts/05_summary.py`)

Writes `results/SUMMARY.md` (using the manuscript names) and collects every table into `results/all_results.xlsx`.
Steps 8 and 9 run later and add their own sections and sheets to these two files.

### Step 6 — Supplementary tables (`scripts/06_supplementary_tables.py`)

Writes `results/Supplementary_Tables_S1_S2.xlsx`: Table S1 (primers from `sample_records/Primers.xlsx`, with the
efficiencies measured here) and Table S2 (per-colony values for the exploratory groups).

### Step 7 — Figures (`scripts/07_figures.py`)

Draws Figures 1 and 2 from the result tables of steps 1, 3 and 4 (no new statistics). Each colony is a point and bars
are group means; effect-size panels show differences with 95 % CI. Fig. 1A marks sex (filled = female, open = male);
Fig. 1B lists p and Benjamini-Hochberg q; Fig. 2B shows the KR repeat run next to the main KR estimate.
Black, bold and * = p < 0.05; grey = p ≥ 0.05. Confidence intervals of the repeat run are in `07_kr_replicate_welch_ci.csv`.

### Step 8 — Copy number v damage (`scripts/08_copy_number_vs_damage.py`)

Question: within a sampling group, do colonies with more estimated mtDNA damage also differ in copy number?
Uses the per-colony values of step 3 (`03_mtdna_per_colony.csv`); no new measurements.

1. **Colonies.** The 27 mtDNA colonies of the main comparisons: nine groups of three (KR, YIS and YKS before and
   during brooding in 2014; coe.3, sp.3 and sp.2 in 2017). Exploratory and archived samples are not used.
   SHC701 (SRP54 failed) is left out of the comparisons that use SRP54 (n = 26).
2. **Measures** (log2, per colony):
   - damage = log2 N0(short COI) − log2 N0(long COI) (= −LS; shown in the figure as lesions per 10 kb);
   - SRP54 per DNA input (the nuclear reference used for copy number);
   - short COI per DNA input (mtDNA quantity without SRP54);
   - SRP54-normalised copy number (short COI − SRP54).
3. **Group effects removed.** Each value is expressed as the difference from its group mean, so the test asks whether
   colonies vary together within groups, not whether groups differ (that is tested in step 3).
4. **Tests.** Null hypothesis: no linear association between the two centred measures. Pearson r with
   df = n − 9 groups − 1; within-group permutation test (values shuffled among colonies of the same group,
   20,000 permutations, seed 1); Spearman ρ (descriptive); r with each group left out in turn.
5. **Result** (`08_copy_number_vs_damage_tests.csv`). Damage was negatively associated with SRP54 per DNA input
   (r = −0.66, df = 16, p = 0.003; permutation p = 0.04; r between −0.78 and −0.45 with any one group left out) and
   not with short COI per DNA input (r = 0.05, df = 17, p = 0.85). SRP54-normalised copy number therefore correlates
   positively with damage (r = 0.60, p = 0.008), but only through SRP54.
6. **Caution.** Damage and short COI per input share the short-COI measurement, so measurement error in short COI
   would push their correlation towards positive values; none was seen. SRP54 and long COI were measured on
   different plates, so their association is not caused by a shared measurement.
7. **Outputs.** `08_copy_number_vs_damage_tests.csv`, `08_copy_number_vs_damage_per_colony.csv` (raw and centred
   values per colony), `FigS1_copy_number_vs_damage.svg/.png`; a "Step 8" section is added to `SUMMARY.md` and two
   sheets to `all_results.xlsx` (running step 8 again replaces them, it does not duplicate them).

### Step 9 — coe.3 in 2014 and 2017 (`scripts/09_coe3_2014_vs_2017.py`)

Question: did coe.3 colonies sampled outside the brooding period give similar damage and copy number in 2014 and 2017,
and how much does each measure vary among coe.3 groups? Uses `03_mtdna_per_colony.csv`; no new measurements.
All colonies of both years were amplified on the same plate for each assay, so their values are directly comparable.

1. **Comparisons.** The three 2017 coe.3 colonies (Yonara, 13 June 2017; about 400 m north of YIS) against
   (A) all nine 2014 coe.3 colonies sampled before brooding (13 July 2014; the reference of Supplementary Table S2) and
   (B) the three YIS colonies sampled before brooding (the nearest 2014 site).
2. **Measures.** Damage (lesions per 10 kb), SRP54-normalised copy number, short COI per DNA input and SRP54 per DNA
   input, as in step 8.
3. **Tests.** Welch's t-test on log2 values; null hypothesis: equal means. 95% CIs from the Welch-Satterthwaite df.
   Not corrected for multiple testing (exploratory).
4. **Result** (`09_coe3_2014_vs_2017_tests.csv`, comparison A). Damage +0.1 lesions per 10 kb (95% CI −4.7 to +4.9;
   p = 0.97); copy number 1.13-fold (p = 0.37); short COI per DNA input 1.27-fold (p = 0.03); SRP54 per DNA input
   1.12-fold (p = 0.44).
5. **Spread among coe.3 groups** (`09_coe3_group_mean_spread.csv`; six 2014 site × period groups and 2017). Group means
   of short COI per DNA input differed by at most 1.38-fold, against 2.14-fold for SRP54-normalised copy number and
   1.96-fold for SRP54 per DNA input.
6. **Caution.** Year is confounded with sampling month (June v July), exact location, sample storage time before the
   2018 qPCR (about 0.8 v 3.7 years) and possibly preservation. This is a consistency check, not a measure of change
   between years.
7. **Outputs.** `09_coe3_2014_vs_2017_tests.csv`, `09_coe3_group_means.csv`, `09_coe3_group_mean_spread.csv`; a
   "Step 9" section in `SUMMARY.md` and three sheets in `all_results.xlsx`.

## 6. How to change something

| To… | Edit |
|---|---|
| use other reference genes for expression | `REFS` in `scripts/01_expression.py` |
| change the outlier rule | the `W["outlier"]` line in `scripts/01_expression.py` |
| exclude or re-include a colony from copy number | `exclude_copy_number` in `data/metadata/samples.csv` |
| change the long-amplicon length | `length_bp` in `qpcr_lib.lesions_per_10kb` |
| change the curve-analysis window | `win_log2` in `qpcr_lib.analyse_curves` (then re-check `02_estimator_validation.csv`) |

## 7. Original instrument files (`eds/` folders)

Each `eds/` folder holds the original StepOnePlus run files (`.eds`), which open in Applied Biosystems StepOne or
QuantStudio Design & Analysis software. The scripts do not read them; they read the `.xls` exports made from the same
runs. Each `.eds` file starts with the same run ID as its exports.

| File | Run | Assays |
|---|---|---|
| `gene_expression/eds/19032018_plate1.eds` | 19 Mar 2018 | cDNA: D2R, CRY1, mtMutS, COI, ACTB, RPL12 |
| `gene_expression/eds/19032018B_plate2.eds` | 19 Mar 2018 | cDNA: SOD, CAT, GPX, HSP70, HSP90, SRP54, ACTB |
| `mtDNA/eds/23032018A_SRP54.eds` | 23 Mar 2018 | DNA: SRP54 (nuclear reference) |
| `mtDNA/eds/23032018B_COIsmall_firstrun.eds` | 23 Mar 2018 | DNA: short COI, first run (replaced; not analysed) |
| `mtDNA/eds/23032018C_COIlarge.eds` | 23 Mar 2018 | DNA: long COI (1,124 bp) |
| `mtDNA/eds/26032018_COIsmall.eds` | 26 Mar 2018 | DNA: short COI, repeat run (analysed) |
| `mtDNA_KR_replicate_20180314/eds/20180314_KR_replicate.eds` | 14 Mar 2018 | DNA: short COI, long COI, SRP54, ACTB (six KR colonies) |

Well names in these files:
- **cDNA plates:** C1-C3 are the colonies sampled before brooding and T1-T3 those sampled during brooding;
  `cDNA_qPCR_platemap.xlsx` and the `expression_label` column of `data/metadata/samples.csv` link them to colony IDs.
- **Main mtDNA plates (23 and 26 March):** wells are named "Sample 1" to "Sample 45". The plate layout (Sheet3) of
  `sample_records/mtDNA_experiment_sample_sheet_20180316.xlsx` gives the colony in each well.
- **14 March replicate run:** wells carry colony IDs (SHC413-SHC415, SHC608, SHC612, SHC614).

## 8. Notes on the data

- Sex was not recorded for the 21 July colonies (SHC593, SHC595, SHC596, SHC608, SHC612, SHC614) or the 2017 panel.
- Reproductive stage was recorded in the field and not confirmed by histology.
- The curve export of the 14 March run lacks wells B4-B9 and F10-F12; instrument Ct values are used for those wells.
- Short-COI melting temperature differed between runs (75.1 °C on 14 March; 77.3 °C on 23 and 26 March 2018);
  ACTB and SRP54 melted identically in both.
- Collection times on 13 July 2014 are approximate.
