# ASD_cortex

Code repository for the manuscript:

**Layered synaptic and developmental-regulatory signatures of ASD risk-gene programs in human cortical transcriptomes**

This repository contains the analysis and figure-generation code used to study how autism spectrum disorder (ASD) risk-gene programs are organized across human cortical development, fetal developmental-state references, adult bulk cortical transcriptomes, and adult single-cell/single-nucleus cortical datasets.

The study asks whether broad ASD risk-gene programs and mid-prenatal-focused subsets represent the same biological signal, or whether they define distinct but related layers of ASD cortical biology.

The main conclusion is a **two-layer model**:

1. **Adult synaptic-dominant layer**  
   A broad ASD risk-gene program, `SFARI_all`, shows the most reproducible adult cortical signal across adult bulk cortical cohorts. This adult signal is largely carried by SynGO-annotated synaptic genes.

2. **Mid-prenatal developmental-regulatory layer**  
   The mid-prenatal-focused program, `midPrenatal_SFARI_top20`, shows sharper fetal developmental localization and stronger chromatin-regulatory context, but does not behave as a stable adult bulk disease anchor.

Together, the analyses support a layered interpretation of ASD cortical molecular convergence: prenatal genetic convergence and adult cortical dysregulation are linked, but the developmental-regulatory and adult synaptic signals are not equivalent.

---

## Repository structure

```text
ASD_cortex/
├── README.md
├── RUN_ORDER.md
├── .gitignore
├── analysis/
│   ├── 01_program_definition/
│   ├── 02_adult_bulk/
│   ├── 03_adult_singlecell/
│   ├── 04_fetal_localization/
│   ├── 05_functional_refinement/
│   ├── 06_specificity/
│   ├── 07_revision_robustness/
│   └── 08_additional_sensitivity/
├── figures/
├── manifests/
│   ├── CODE_FILE_MANIFEST.tsv
│   └── ZENODO_INPUTS_TEMPLATE.tsv
└── env/
    └── sessionInfo_generated.txt
```

Large processed data objects, generated results, runtime logs, and submission scripts are not included in this repository. Processed input files, harmonized intermediate tables, source tables, and supplementary tables supporting the analyses are available from Zenodo:

**https://doi.org/10.5281/zenodo.20046256**

---

## What is included

This repository includes analysis code only.

Included:

- R scripts for program definition, adult bulk scoring, fetal localization, adult single-cell localization, deconvolution, enrichment analysis, robustness analysis, and final sensitivity analyses.
- Figure-generation scripts.
- A code manifest with file sizes and MD5 checksums.
- A Zenodo input template describing the expected processed intermediate files.
- R session information, if available.

Excluded:

- Large processed input objects.
- Runtime logs.
- Generated result tables.
- Generated figures.
- Raw public datasets.

The processed data package is deposited in Zenodo and should be used together with this code repository.

---

## Data availability

Processed input files, harmonized intermediate tables, source tables, and supplementary tables supporting the analyses are available from Zenodo:

**https://doi.org/10.5281/zenodo.20046256**

The Zenodo archive contains processed and harmonized files needed for practical reproduction of the manuscript analyses. Raw public datasets are not redistributed in the archive and should be obtained from their original repositories.

Expected Zenodo input organization is summarized in:

```text
manifests/ZENODO_INPUTS_TEMPLATE.tsv
```

---

## Main data resources

The analyses use publicly available and de-identified human datasets, including:

- **BrainSpan** human brain developmental transcriptomic reference;
- **Eze 2021** early human brain single-cell developmental-state reference;
- **Nowakowski-UCSC** fetal cortex reference;
- **GSE162170** fetal RNA/multiome reference;
- **GSE102741** adult bulk cortical transcriptomic cohort;
- **GSE64018** adult bulk cortical transcriptomic cohort;
- **Gandal2022** postmortem adult cortical transcriptomic resource;
- **Velmeshev** adult cortical single-nucleus dataset;
- **PsychENCODE 2024** adult cortical ASD single-cell/single-nucleus resource;
- **SynGO** synaptic ontology/resource;
- **Reactome / MSigDB** gene-set resources.


---

## Analysis overview

### 1. ASD risk-gene program definition

The study begins with a standardized primary SFARI source pool. Four predefined ASD risk-gene programs are used throughout the manuscript:

- `SFARI_all`
- `midPrenatal_SFARI_top20`
- `midPrenatal_SFARI_top10`
- `midPrenatal_SFARI_top05`

The broad program `SFARI_all` is defined as the intersection between the standardized primary SFARI source pool and the BrainSpan-derived core gene universe with non-missing mid-prenatal ranking values.

The mid-prenatal programs are derived using a BrainSpan cortical developmental concentration metric. In brief, genes are ranked according to how strongly their mid-prenatal cortical expression is elevated relative to their own across-stage cortical background.

Relevant scripts:

```text
analysis/01_program_definition/
```

Key expected outputs:

- standardized SFARI source-pool gene list;
- predefined program membership table;
- BrainSpan developmental ranking;
- program-size and overlap summaries.

---

### 2. Adult bulk cortical program scoring

Adult bulk cortical analyses are performed in:

- GSE102741
- GSE64018
- Gandal2022

For each cohort, predefined programs are scored using mapped member genes. The primary model compares ASD versus control program scores using diagnosis-based linear models. Cross-cohort summaries are used to evaluate directionality and consistency.

Relevant scripts:

```text
analysis/02_adult_bulk/
```

Key expected outputs:

- sample-level program scores;
- cohort-level ASD-control model summaries;
- cross-cohort direction summaries;
- random-effects meta-analysis sensitivity summaries.

---

### 3. Adult single-cell and single-nucleus localization

Adult single-cell/single-nucleus analyses are performed using:

- PsychENCODE
- Velmeshev

The main goal is to identify the adult cortical cellular context of the ASD risk-gene programs. Program localization is summarized across broad cell classes, especially excitatory and inhibitory neuronal compartments.

Relevant scripts:

```text
analysis/03_adult_singlecell/
```

Key expected outputs:

- broad-class program localization summaries;
- UCell/AUCell sensitivity summaries;
- donor-level mean-score and pseudobulk summaries;
- donor-aware case-control model summaries.

---

### 4. Fetal developmental-state localization

Fetal developmental-state localization analyses test whether the ASD risk-gene programs localize to specific developmental states, including progenitor-to-neurogenic contexts.

The primary fetal reference is the Eze 2021 early human brain single-cell dataset. Independent fetal cortical references include Nowakowski-UCSC and GSE162170 RNA/multiome resources.

Relevant scripts:

```text
analysis/04_fetal_localization/
```

Key expected outputs:

- fetal cluster-stage program summaries;
- broad fetal class localization summaries;
- top fetal developmental-state tables;
- independent fetal reference summaries.

---

### 5. Functional refinement

Functional refinement analyses are used to distinguish synaptic and regulatory features of the ASD risk-gene programs.

The analyses include:

- SynGO enrichment;
- Reactome enrichment;
- transcription-factor target enrichment;
- SynGO-annotated synaptic versus non-synaptic decomposition.

Relevant scripts:

```text
analysis/05_functional_refinement/
```

Key expected outputs:

- SynGO enrichment tables;
- Reactome enrichment tables;
- transcription-factor target enrichment tables;
- synaptic and non-synaptic decomposition summaries.

---

### 6. Specificity analyses

Specificity analyses test whether observed enrichment patterns exceed matched-random expectations.

The analyses include:

- size-matched random controls;
- expression-matched random controls;
- empirical P-value summaries;
- standardized null-deviation summaries.

Relevant scripts:

```text
analysis/06_specificity/
```

Key expected outputs:

- matched-random empirical P-value summaries;
- matched-random z-score summaries;
- integrated specificity summaries.

---

### 7. Robustness analyses

Additional robustness analyses support the main interpretation and address cohort-specific behavior.

These include:

- alternative adult bulk scoring methods;
- random-effects adult-bulk meta-analysis;
- gene-set-definition sensitivity;
- dataset-specific gene-harmonization summaries;
- GSE64018 divergence and sample-influence analysis;
- synaptic leave-one-gene-out analysis;
- top-expression-dropout analysis;
- best-available covariate-adjusted model summaries;
- CAMERA and mROAST inter-gene-correlation-aware sensitivity analyses;
- alternative deconvolution analyses.

Relevant scripts:

```text
analysis/07_revision_robustness/
```


---

### 8. Final additional sensitivity analyses

Two final sensitivity analyses were added to strengthen the final manuscript version.

#### 8.1 Detection-filtered adult bulk rescoring

Genes were retained if detected in at least 50% of samples in each cohort. Mapped-all and detected-only program scores were compared for:

- `SFARI_all`
- `SFARI_all_synaptic`
- `midPrenatal_SFARI_top20`

The analysis showed that detection filtering preserved the direction of `SFARI_all` and `SFARI_all_synaptic` effects across adult bulk cohorts, while `midPrenatal_SFARI_top20` remained cohort-dependent.

Relevant scripts:

```text
analysis/08_additional_sensitivity/Step60B_detection_filtered_rescoring_v3.R
analysis/08_additional_sensitivity/Step60B_add_Gandal_and_merge_v4.R
```

Expected final table in the Zenodo archive:

```text
32_manuscript_compact_three_cohort_detection_summary_v4.tsv
```

#### 8.2 BrainSpan prenatal-boundary sensitivity

SFARI source-pool genes were re-ranked using alternative BrainSpan prenatal boundary definitions:

- early-prenatal;
- mid-prenatal;
- late-prenatal;
- early-to-mid prenatal;
- mid-to-late prenatal;
- all-prenatal;
- mid-minus-adjacent prenatal specificity.

The primary mid-prenatal ranking recovered the predefined `midPrenatal_SFARI_top20` set. Adjacent or broader prenatal definitions produced partly overlapping gene sets and did not convert the mid-prenatal program into a stable adult bulk anchor.

Relevant script:

```text
analysis/08_additional_sensitivity/Step60A_BrainSpan_boundary_sensitivity_v1.R
```

Expected final table in the Zenodo archive:

```text
09_manuscript_compact_boundary_sensitivity_summary.tsv
```

---

## Figure generation

Figure-generation scripts are stored in:

```text
figures/
```

The main figures correspond to:

- **Figure 1:** integrated program prioritization across developmental concentration and adult bulk support;
- **Figure 2:** adult cortical support and neuronal-context localization;
- **Figure 3:** bulk composition remodeling and composition-adjusted interpretation;
- **Figure 4:** fetal developmental-state localization;
- **Figure 5:** SynGO and Reactome functional refinement;
- **Figure 6:** matched-random specificity analyses.

Generated figures are not stored in this code repository. Source tables needed for figure reproduction are included in the Zenodo archive.

---

## How to reproduce the analyses

This repository is intended to be used together with the Zenodo processed-data archive.

Recommended steps:

1. Clone this repository.

```bash
git clone https://github.com/aojielian/ASD_cortex.git
cd ASD_cortex
```

2. Download the processed Zenodo data package:

```text
https://doi.org/10.5281/zenodo.20046256
```

3. Arrange processed input files according to:

```text
manifests/ZENODO_INPUTS_TEMPLATE.tsv
```

4. Check or edit file paths in each script.

5. Run scripts in the order described in:

```text
RUN_ORDER.md
```

The scripts were written as modular R scripts rather than a single monolithic workflow. This makes it easier to inspect, rerun, and validate each analysis module independently.

---

## Software requirements

The analyses were performed primarily in R.

Main R packages include:

- `data.table`
- `Matrix`
- `ggplot2`
- `limma`
- `nnls`
- `openxlsx`
- `AnnotationDbi`
- `org.Hs.eg.db`
- `Seurat`
- `UCell`
- `GSVA`

Additional packages may be required depending on the specific analysis module.

R session information is provided in:

```text
env/sessionInfo_generated.txt
```

---

## Reproducibility notes

Large intermediate objects are intentionally excluded from GitHub. Processed data are hosted through Zenodo. Slurm submission scripts and runtime logs are not included. Generated figures and result tables are also not included in the code repository.

The code was organized for transparency and traceability. It is not intended to be a fully containerized one-command workflow.

---

## Citation

If you use this repository, please cite the associated manuscript and data archive.

Manuscript:

**Layered synaptic and developmental-regulatory signatures of ASD risk-gene programs in human cortical transcriptomes**

Processed data archive:

Lian A. Processed data and source tables for “Layered synaptic and developmental-regulatory signatures of ASD risk-gene programs in human cortical transcriptomes”. Zenodo. 2026. https://doi.org/10.5281/zenodo.20046256

Citation details for the manuscript will be added after publication.

---

## Contact

Aojie Lian, PhD  
NHC Key Laboratory of Birth Defect for Research and Prevention  
Hunan Provincial Maternal and Child Health Care Hospital  
Changsha, Hunan 410028, China  
Email: aojielian@gmail.com
