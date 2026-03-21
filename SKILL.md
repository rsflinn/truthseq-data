---
name: truthseq-validate
description: Validate gene regulatory predictions against real experimental data from 11,000 human gene knockdowns. Clone the repo, set up the Replogle Perturb-seq reference data, then test user-supplied or example claims across up to three tiers of evidence.
allowed-tools: Bash(git *), Bash(python3 *), Bash(pip3 *), Bash(curl *), Bash(cat *), Bash(ls *), Bash(head *), Bash(wc *)
---

# TruthSeq: Validate Gene Regulatory Claims Against Lab Data

## Overview

TruthSeq checks whether computational predictions about gene regulation hold up against real experimental data. Given a CSV of claims ("Gene X regulates Gene Y in direction Z"), it tests each claim against:

- **Tier 1**: The Replogle Perturb-seq atlas (~11,000 single-gene CRISPR knockdowns in human K562 cells)
- **Tier 2** (optional): Disease tissue gene expression data from public repositories
- **Tier 3** (optional): Genetic association scores from Open Targets

Each claim receives a grade: VALIDATED, PARTIALLY_SUPPORTED, WEAK, CONTRADICTED, or UNTESTABLE.

## Prerequisites

- Python 3.8+
- ~500 MB free disk space
- Internet connection (for Figshare download and optional Open Targets queries)

## Step 1: Clone the repository

```bash
git clone https://github.com/rsflinn/truthseq.git
cd truthseq
```

## Step 2: Install dependencies

```bash
pip3 install scanpy anndata pandas pyarrow numpy scipy requests
```

## Step 3: Download and process the reference dataset

```bash
python3 setup.py --skip-gene-map
```

This downloads the Replogle Perturb-seq pseudo-bulk h5ad (~357 MB) from Figshare and processes it into two parquet files:
- `replogle_knockdown_effects.parquet` — ~37.7 million gene-gene pairs with Z-scores
- `replogle_knockdown_stats.parquet` — per-knockdown distribution statistics for percentile calculations

Verify setup succeeded:

```bash
python3 setup.py --status
```

Expected: both parquet files show `[OK]`.

## Step 4: Run example validation (Tier 1 only)

```bash
python3 truthseq_validate.py \
    --claims example_claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --output example_results
```

### Expected output

The example file contains 11 claims designed to produce all five grade types:

| Claim type | Expected grade | Count |
|-----------|---------------|-------|
| Known biology (correct direction) | PARTIALLY_SUPPORTED | 5 |
| Wrong direction controls | CONTRADICTED | 2 |
| Weak/absent signal controls | WEAK | 3 |
| Gene not in dataset | UNTESTABLE | 1 |

Note: PARTIALLY_SUPPORTED (not VALIDATED) is the highest possible grade without Tier 2 disease data.

### Verify results

```bash
cat example_results/truthseq_results.csv
```

Check that:
- SLC30A1 → MT2A scores PARTIALLY_SUPPORTED with Z-score ~76
- GATA1 → TYROBP scores PARTIALLY_SUPPORTED with Z-score ~30
- SLC30A1 → MT2A (UP) scores CONTRADICTED (same pair, wrong direction)
- FOXP2 → CNTNAP2 scores UNTESTABLE (FOXP2 not in knockdown dataset)

## Step 5: Run with disease context (Tier 2)

To reach VALIDATED, claims need both Tier 1 perturbation evidence AND Tier 2 disease tissue evidence. Supply a disease expression file:

```bash
python3 truthseq_validate.py \
    --claims your_claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --disease-expr your_disease_de_results.tsv \
    --output results_with_disease
```

Disease expression files should have columns: `gene`, `log2fc`, `padj` (and optionally `cell_type`). See `format_spec.md` for details.

Alternatively, search for publicly available datasets:

```bash
python3 dataset_search.py --query "autism brain RNA-seq" --verbose
```

Or let TruthSeq search its built-in registry automatically:

```bash
python3 truthseq_validate.py \
    --claims your_claims.csv \
    --disease "breast cancer" \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet
```

## Step 6: Create your own claims file

Create a CSV with these columns:

```csv
upstream_gene,downstream_gene,predicted_direction,cell_type_context,source
MYT1L,MEF2C,UP,neuron,GRN inference
TP53,CDKN1A,UP,,literature
GATA1,HBB,UP,,my_analysis
```

- `upstream_gene`: the predicted regulator
- `downstream_gene`: the predicted target
- `predicted_direction`: UP means the regulator activates the target when active; DOWN means it represses it
- `cell_type_context`: optional, for annotation
- `source`: optional, for your tracking

Run validation:

```bash
python3 truthseq_validate.py \
    --claims my_claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --output my_results
```

## Step 7: Interpret the report

The output directory contains:
- `truthseq_results.csv` — full evidence table with Z-scores, percentiles, disease expression, and grades
- `truthseq_summary.md` — human-readable report
- `truthseq_heatmap.png` — visual summary of confidence grades

### Grade definitions

- **VALIDATED**: Tier 1 knockdown confirms direction (top 10% effect) AND Tier 2 disease tissue shows the target gene is significantly dysregulated
- **PARTIALLY_SUPPORTED**: Tier 1 confirms direction but disease data is missing or non-significant
- **WEAK**: The regulator was tested but the target gene didn't respond notably
- **CONTRADICTED**: Tier 1 shows the opposite direction from the prediction
- **UNTESTABLE**: The regulator gene isn't in the knockdown dataset

## Caveats

- Tier 1 data is from K562 cells (blood-derived). Gene regulation varies by cell type. A WEAK grade means "not detectable in this cell type," not "definitely wrong."
- The dataset contains ~7,600 unique knockdowns covering ~8,200 target genes. Not all human genes are represented.
- Z-scores are normalized within each knockdown. Highly pleiotropic genes (affecting thousands of targets) will have individual effects diluted.
- VALIDATED requires both perturbation and disease evidence. Without Tier 2 data, the maximum grade is PARTIALLY_SUPPORTED.
