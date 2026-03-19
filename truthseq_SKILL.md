---
name: truthseq
description: >
  TruthSeq validates computational gene-gene regulatory claims against real
  perturbation experiments. Given a set of predicted regulatory relationships
  (e.g., "MYT1L knockdown reduces MEF2C expression"), TruthSeq checks whether
  those predictions hold up in actual CRISPR perturbation data, disease tissue
  expression, and genetic association databases. Outputs a confidence score
  for each claim. Use when you want to reality-check in silico genomics
  findings against in vivo or ex vivo experimental evidence.
allowed-tools: Bash(python3 *), Bash(pip *), Bash(curl *), Bash(wget *), WebFetch
---

# TruthSeq — From In Silico to In Vivo

Validate computational gene regulatory claims against real experimental data.

## Overview

TruthSeq takes a table of gene-gene regulatory predictions from any
computational analysis (network inference, pathway enrichment, co-expression,
etc.) and checks each prediction against:

1. **Direct perturbation data** from genome-wide CRISPR screens (Replogle et al. 2022)
2. **Disease tissue expression** from postmortem brain single-cell RNA-seq (PsychENCODE)
3. **Genetic associations** from the Open Targets Platform (live API)

Each claim receives a confidence grade: VALIDATED, PARTIALLY SUPPORTED,
CONTRADICTED, WEAK, or UNTESTABLE — with full evidence and base-rate context.

## Prerequisites

Python 3.9+ with pip. The skill will install its own dependencies.

## Step 1: Set up environment

```bash
pip install pandas pyarrow requests scipy matplotlib seaborn --quiet
```

## Step 2: Prepare input claims

Create a CSV file called `claims.csv` with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| upstream_gene | Yes | Gene symbol of the claimed regulator (e.g., MYT1L) |
| downstream_gene | Yes | Gene symbol of the claimed target (e.g., MEF2C) |
| predicted_direction | Yes | DOWN or UP — predicted effect on downstream when upstream is knocked down |
| cell_type_context | No | Tissue/cell type context (e.g., "neuron", "excitatory_neuron") |
| source | No | Where this claim came from (e.g., "GRN inference", "pathway analysis") |

**Example claims.csv** (from an autism gene regulatory network hypothesis):
```csv
upstream_gene,downstream_gene,predicted_direction,cell_type_context,source
MYT1L,MEF2C,DOWN,neuron,GRN inference
TCF4,MEF2C,DOWN,neuron,GRN inference
MYT1L,SCN2A,DOWN,excitatory_neuron,co-expression
EP300,MYT1L,DOWN,neuron,chromatin regulation
FOXP1,TCF4,DOWN,neuron,TF binding prediction
MEF2C,KCNA2,DOWN,neuron,target gene analysis
MEF2C,GRIN2B,DOWN,excitatory_neuron,target gene analysis
MECP2,CDKL5,DOWN,neuron,literature
```

If you do not have a claims file, the skill will generate the example above
as a demonstration.

## Step 3: Download reference datasets

Download the pre-processed perturbation and expression reference data:

```bash
# Replogle Perturb-seq atlas — pre-processed knockdown effects
# ~9,867 genes knocked down in K562 cells, effect on all other genes
wget -q -O replogle_effects.parquet \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/replogle_knockdown_effects.parquet"

# PsychENCODE ASD differential expression by cell type
wget -q -O psychencode_de.parquet \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/psychencode_asd_de.parquet"

# Gene symbol to Ensembl ID mapping
wget -q -O gene_map.tsv \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/gene_id_mapping.tsv"
```

> **Note:** If the hosted data files are unavailable, the skill falls back to
> querying the Open Targets API directly and reports which validation tiers
> could not be checked due to missing reference data.

## Step 4: Download TruthSeq validation script

```bash
wget -q -O truthseq_validate.py \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/truthseq_validate.py"
```

## Step 5: Run TruthSeq validation

```bash
python3 truthseq_validate.py \
  --claims claims.csv \
  --replogle replogle_effects.parquet \
  --psychencode psychencode_de.parquet \
  --gene-map gene_map.tsv \
  --output truthseq_report
```

This script performs the following:

### 5a. Gene validation and ID mapping
- Validates all gene symbols against the reference mapping
- Flags unrecognized genes
- Maps symbols to Ensembl IDs for API queries

### 5b. Perturbation lookup (Tier 1 — direct experimental evidence)
For each upstream->downstream claim:
- Checks if upstream_gene was knocked down in the Replogle Perturb-seq atlas
- Retrieves the observed log2 fold change of downstream_gene after knockdown
- Computes a null distribution: the log2FC of 500 random genes after the same
  knockdown, to establish a base rate
- Calculates the percentile rank of the downstream gene's response vs null
- Compares observed direction to predicted direction

### 5c. Disease tissue expression check (Tier 2 — observational)
For each downstream_gene:
- Looks up cell-type-specific differential expression in ASD vs control brain
  (PsychENCODE single-nucleus RNA-seq)
- Reports whether the gene is actually dysregulated in disease tissue
- If cell_type_context is specified, checks that specific cell type
- Flags when data exists only in a different cell type than claimed

### 5d. Open Targets query (Tier 3 — genetic association)
For each gene:
- Queries the Open Targets Platform GraphQL API for disease associations
- Reports known genetic links to autism spectrum disorder or epilepsy
- This provides population-level context, not causal validation

### 5e. Confidence scoring
Each claim receives a composite grade:

| Grade | Criteria |
|-------|----------|
| **VALIDATED** | Perturbation data confirms predicted direction AND effect ranks >90th percentile vs null AND disease tissue is consistent |
| **PARTIALLY SUPPORTED** | Effect is real but modest (50-90th percentile), OR direction matches only in a different cell type |
| **CONTRADICTED** | Perturbation data shows significant effect in the OPPOSITE direction |
| **WEAK** | Perturbation data exists but downstream gene is not more affected than random genes |
| **UNTESTABLE** | No perturbation data for the upstream gene in any available dataset |
| **CELL-TYPE CAVEAT** | Evidence exists but only in a cell type different from the claimed context |

## Step 6: Review outputs

The skill produces three output files in the `truthseq_report/` directory:

### truthseq_results.csv
Full evidence table with one row per claim:
```
upstream_gene, downstream_gene, predicted_direction, perturb_log2fc,
perturb_pvalue, perturb_percentile, perturb_direction_match,
perturb_cell_type, psychencode_log2fc, psychencode_padj,
psychencode_cell_type, ot_association_score, confidence_grade, evidence_summary
```

### truthseq_summary.md
Human-readable report:
- Overall scorecard (X validated, Y contradicted, Z untestable)
- Per-claim evidence narratives
- Base-rate context (what fraction of random gene pairs would score this well?)
- Cell-type coverage warnings
- Recommendations for experimental follow-up

### truthseq_heatmap.png
Visualization showing claims x evidence sources, color-coded by confidence grade.

## Expected output (example)

For the demonstration ORC gene claims, a well-functioning run produces:

```
============================================================
TruthSeq Validation Complete
============================================================
  VALIDATED                 0
  PARTIALLY_SUPPORTED       7
  CELL_TYPE_CAVEAT          0
  WEAK                      1
  CONTRADICTED              0
  UNTESTABLE                0
  TOTAL                     8
```

Key findings from the demo:
- 7 claims are PARTIALLY_SUPPORTED: no perturbation effect in K562 cells
  (expected — these are neuron-specific regulatory relationships being tested
  in a leukemia cell line), but downstream genes ARE dysregulated in ASD
  postmortem brain tissue.
- 1 claim is WEAK (MECP2 -> CDKL5): the tool correctly identifies that the
  predicted regulatory direction is wrong — CDKL5 phosphorylates MECP2,
  not the reverse.
- The cell-type mismatch warning highlights that neuronal perturbation
  datasets (Wu 2022, Tian 2019) would be needed to validate neuron-specific
  claims.

## Fallback mode (no pre-processed data)

If reference data files cannot be downloaded, the skill falls back to:
1. Query Open Targets API only (Tier 3)
2. Report which tiers are missing and what data would be needed
3. Still produces a partial confidence report with available evidence
4. Explicitly states: "Tier 1 (perturbation) and Tier 2 (disease tissue)
   validation could not be performed. Results reflect genetic association
   evidence only."

## Reproducibility notes

- All reference data is derived from published, publicly available sources
- Replogle data: GEO accession GSE132080 (Replogle et al., Cell 2022)
- PsychENCODE data: PsychENCODE Consortium portal (Gandal et al., Science 2018;
  updated 2022 release)
- Open Targets: api.platform.opentargets.org (queried live at runtime)
- Pre-processing scripts are included in the truthseq-data GitHub repository
- Random seed is fixed (42) for null distribution reproducibility
