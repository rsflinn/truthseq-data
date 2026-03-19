# TruthSeq

Validate computational gene regulatory claims against real perturbation experiments.

## What it does

TruthSeq takes a table of gene-gene regulatory predictions — the kind produced by network inference, pathway enrichment, co-expression analysis or any other computational method — and checks each one against experimental data from public databases. It tells you which claims survive contact with actual biology and which don't.

For each claim, TruthSeq reports a confidence grade (VALIDATED, PARTIALLY SUPPORTED, CONTRADICTED, WEAK, or UNTESTABLE) with full evidence and base-rate context showing how the claimed gene pair compares to random pairs from the same dataset.

## Why it exists

Computational genomics tools make it easy to generate thousands of gene regulatory predictions that look statistically significant but live entirely in annotation space. Most are never checked against what actually happens when you knock a gene down in a real cell.

TruthSeq was built during a citizen-science project investigating autism genetics, where over 90 computational analyses produced a causal hypothesis supported by multiple standard network metrics. When tested against actual perturbation data, the claim did not survive (p = 0.31). TruthSeq is designed to catch that kind of failure early.

## Validation tiers

**Tier 1 — Direct perturbation (gold standard).** The Replogle et al. (2022) Perturb-seq atlas provides genome-wide CRISPRi knockdown data for ~9,867 genes in K562 cells. For each claim, TruthSeq looks up whether the upstream gene was knocked down and what happened to the downstream gene's expression. It computes a null distribution from 500 random genes' responses to the same knockdown, so the output answers not just "did the downstream gene change?" but "did it change more than a random gene would?" — the base-rate question most analyses skip.

**Tier 2 — Disease tissue expression (observational).** PsychENCODE single-nucleus RNA-seq from ASD vs. control postmortem brain provides cell-type-specific differential expression. TruthSeq checks whether each downstream gene is actually dysregulated in the relevant disease context and flags cell-type mismatches.

**Tier 3 — Genetic association (population-level).** The Open Targets Platform API provides gene-disease association scores for population-level context.

## Confidence grades

| Grade | Meaning |
|---|---|
| **VALIDATED** | Perturbation confirms predicted direction (>90th percentile vs. null) AND disease tissue is consistent |
| **PARTIALLY SUPPORTED** | Effect is real but modest (50th-90th percentile), or direction matches only in a different cell type |
| **CONTRADICTED** | Perturbation data shows significant effect in the opposite direction |
| **WEAK** | Perturbation data exists but downstream gene is not more affected than random genes |
| **UNTESTABLE** | No perturbation data available for the upstream gene |
| **CELL-TYPE CAVEAT** | Evidence exists but only in a different cell type than claimed |

## Quick start

### Requirements

Python 3.9+ with pip.

### Install dependencies

```bash
pip install pandas pyarrow requests scipy matplotlib seaborn
```

### Download reference data

```bash
# Replogle Perturb-seq atlas — pre-processed knockdown effects (~9,867 genes)
wget -q -O replogle_knockdown_effects.parquet \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/replogle_knockdown_effects.parquet"

# PsychENCODE ASD differential expression by cell type
wget -q -O psychencode_asd_de.parquet \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/psychencode_asd_de.parquet"

# Gene symbol to Ensembl ID mapping
wget -q -O gene_id_mapping.tsv \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/gene_id_mapping.tsv"

# Validation script
wget -q -O truthseq_validate.py \
  "https://github.com/rsflinn/truthseq-data/releases/download/v1.0/truthseq_validate.py"
```

### Prepare your claims

Create a CSV file with your regulatory predictions:

```csv
upstream_gene,downstream_gene,predicted_direction,cell_type_context,source
MYT1L,MEF2C,DOWN,neuron,GRN inference
TCF4,MEF2C,DOWN,neuron,GRN inference
EP300,MYT1L,DOWN,neuron,chromatin regulation
MEF2C,KCNA2,DOWN,neuron,target gene analysis
MECP2,CDKL5,DOWN,neuron,literature
```

| Column | Required | Description |
|---|---|---|
| upstream_gene | Yes | Gene symbol of the claimed regulator |
| downstream_gene | Yes | Gene symbol of the claimed target |
| predicted_direction | Yes | DOWN or UP — predicted effect on downstream when upstream is knocked down |
| cell_type_context | No | Tissue or cell type context (e.g., "neuron", "excitatory_neuron") |
| source | No | Where this claim came from |

If you run the script without a claims file, it uses a built-in demo set.

### Run

```bash
python3 truthseq_validate.py \
  --claims claims.csv \
  --replogle replogle_knockdown_effects.parquet \
  --psychencode psychencode_asd_de.parquet \
  --gene-map gene_id_mapping.tsv \
  --output truthseq_report
```

Options:
- `--skip-ot` skips Open Targets API queries (faster, offline)
- `--skip-base-rate` skips the 1,000-iteration null simulation (faster)

### Output

The script produces three files in the output directory:

- `truthseq_results.csv` — Full evidence table with one row per claim, all scores and evidence narratives
- `truthseq_summary.md` — Human-readable report with overall scorecard, per-claim evidence and base-rate comparison
- `truthseq_heatmap.png` — Visual summary color-coded by confidence grade

## Example output

Running the demo claims (8 regulatory predictions from an autism gene network hypothesis):

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

Seven claims are PARTIALLY SUPPORTED: no perturbation effect in K562 cells (expected — these are neuron-specific relationships tested in a leukemia cell line), but the downstream genes are dysregulated in ASD postmortem brain. One claim is WEAK (MECP2 → CDKL5): the tool correctly flags that the predicted regulatory direction is unsupported — CDKL5 phosphorylates MECP2, not the reverse.

## Fallback mode

If the reference data files can't be downloaded, TruthSeq falls back to querying the Open Targets API only (Tier 3) and reports which tiers couldn't be checked. The output explicitly states what data would be needed for full validation.

## Limitations

**Cell-type mismatch** is the primary limitation. The Replogle atlas uses K562 cells (a leukemia line), not neurons. Gene regulatory relationships differ across cell types, and TruthSeq flags this mismatch rather than hiding it. As neuronal Perturb-seq datasets grow, the reference data can be updated.

**Pairwise only.** TruthSeq validates individual gene-gene claims, not network-level properties like cascade amplification or multi-step propagation. Extending to multi-hop validation is a future direction.

**Perturbation ≠ direct regulation.** A knockdown effect shows that gene X's expression affects gene Y, but the mechanism could be indirect.

## AI agent integration

TruthSeq is designed to be executed by AI coding agents. A SKILL.md file is included that provides step-by-step instructions any agent can follow to download the data, run the validation and interpret results. See [`truthseq_SKILL.md`](truthseq_SKILL.md) in this repository.

## Data sources

All reference data is derived from published, publicly available sources:

- **Replogle Perturb-seq atlas**: Replogle JM, et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell*. 2022;185(14):2559-2575. GEO: GSE132080
- **PsychENCODE**: Gandal MJ, et al. Transcriptome-wide isoform-level dysregulation in ASD, schizophrenia, and bipolar disorder. *Science*. 2018;362(6420):eaat8127
- **Open Targets Platform**: https://platform.opentargets.org (queried live at runtime)

## License

MIT
