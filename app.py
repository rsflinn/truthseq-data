"""
TruthSeq Web Interface
======================
A Streamlit app that lets anyone validate gene regulatory claims
without installing Python or running terminal commands.

Run locally:
    streamlit run app.py

Deploy on Streamlit Community Cloud:
    1. Push this file to your GitHub repo
    2. Go to share.streamlit.io
    3. Point it at your repo and app.py
    4. The app will download data on first run (~500 MB)

Requirements:
    pip install streamlit pandas pyarrow numpy scipy
"""

import os
import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats as sp_stats
from io import StringIO

# ============================================================
# Configuration
# ============================================================

PERCENTILE_VALIDATED = 90
PERCENTILE_PARTIAL = 50

GRADE_COLORS = {
    'VALIDATED': '#2ecc71',
    'PARTIALLY_SUPPORTED': '#f1c40f',
    'WEAK': '#bdc3c7',
    'CONTRADICTED': '#e74c3c',
    'UNTESTABLE': '#95a5a6',
}

GRADE_DESCRIPTIONS = {
    'VALIDATED': 'Lab data confirms the direction and strength of this regulatory relationship.',
    'PARTIALLY_SUPPORTED': 'Direction matches but effect size is modest, or only partial evidence available.',
    'WEAK': 'The upstream gene was tested but the downstream gene did not respond strongly. '
            'This does not mean the relationship is wrong — it may be cell-type-specific.',
    'CONTRADICTED': 'Lab data shows the opposite direction from what was predicted.',
    'UNTESTABLE': 'The upstream gene was not in the knockdown dataset.',
}

# ============================================================
# Data Loading (cached)
# ============================================================

@st.cache_data(show_spinner="Loading knockdown data...")
def load_replogle_data(effects_path, stats_path=None):
    """Load the Replogle Perturb-seq data."""
    effects_df = pd.read_parquet(effects_path)
    stats_df = None
    if stats_path and os.path.exists(stats_path):
        stats_df = pd.read_parquet(stats_path)
    return effects_df, stats_df


def find_data_files():
    """Look for parquet files in the current directory and common locations."""
    search_dirs = ['.', 'data', os.path.expanduser('~')]
    found = {}

    for search_dir in search_dirs:
        for fname in os.listdir(search_dir) if os.path.isdir(search_dir) else []:
            fpath = os.path.join(search_dir, fname)
            if fname.endswith('_knockdown_effects.parquet'):
                ct = 'K562' if 'rpe1' not in fname.lower() else 'RPE1'
                found[f'{ct}_effects'] = fpath
            elif fname.endswith('_knockdown_stats.parquet'):
                ct = 'K562' if 'rpe1' not in fname.lower() else 'RPE1'
                found[f'{ct}_stats'] = fpath

    return found


# ============================================================
# Validation Logic (simplified from truthseq_validate.py)
# ============================================================

def compute_percentile_from_distribution(z_score, kd_data):
    """Compute percentile rank of a Z-score within a knockdown's distribution."""
    abs_z = abs(z_score)
    all_abs_z = kd_data['z_score'].abs().values
    return float(sp_stats.percentileofscore(all_abs_z, abs_z))


def validate_single_claim(upstream, downstream, predicted_dir, effects_df, stats_df=None):
    """Validate a single gene regulatory claim."""
    available_kd = set(effects_df['knocked_down_gene'].unique())
    stats_kd = set(stats_df['knocked_down_gene'].unique()) if stats_df is not None else set()
    all_kd = available_kd | stats_kd

    if upstream not in all_kd:
        return {
            'grade': 'UNTESTABLE',
            'z_score': None,
            'percentile': None,
            'direction_match': None,
            'cell_line': None,
            'note': f'{upstream} was not in the knockdown dataset.',
        }

    # Look up the knockdown effect
    kd_data = effects_df[effects_df['knocked_down_gene'] == upstream]
    target_hit = kd_data[kd_data['affected_gene'] == downstream]

    if len(target_hit) > 0:
        # Per-cell-type results
        cell_types = target_hit['cell_line'].unique() if 'cell_line' in target_hit.columns else ['K562']
        per_ct = []

        for ct in cell_types:
            ct_hits = target_hit[target_hit['cell_line'] == ct] if 'cell_line' in target_hit.columns else target_hit
            z = float(ct_hits.iloc[0]['z_score'])
            obs_dir = "UP" if z < 0 else "DOWN"
            dir_match = (obs_dir == predicted_dir)

            ct_kd = kd_data
            if 'cell_line' in kd_data.columns:
                ct_kd = kd_data[kd_data['cell_line'] == ct]

            pct = compute_percentile_from_distribution(z, ct_kd) if len(ct_kd) > 0 else 50.0

            per_ct.append({
                'cell_line': ct,
                'z_score': round(z, 4),
                'percentile': round(pct, 1),
                'direction_match': dir_match,
                'observed_direction': obs_dir,
            })

        # Use best result
        best = max(per_ct, key=lambda x: x['percentile'])
        z = best['z_score']
        pct = best['percentile']
        dir_match = best['direction_match']
        cell_line = best['cell_line']

        # Assign grade
        if dir_match and pct >= PERCENTILE_VALIDATED:
            grade = 'VALIDATED'
        elif dir_match and pct >= PERCENTILE_PARTIAL:
            grade = 'PARTIALLY_SUPPORTED'
        elif not dir_match and pct >= PERCENTILE_PARTIAL:
            grade = 'CONTRADICTED'
        else:
            grade = 'WEAK'

        note = f"Z={z:.2f}, {pct:.0f}th percentile in {cell_line}"
        if len(per_ct) > 1:
            ct_notes = [f"{r['cell_line']}: Z={r['z_score']:.2f}, {r['percentile']:.0f}th pctl"
                        for r in per_ct]
            note += f". All cell types: {'; '.join(ct_notes)}"

        return {
            'grade': grade,
            'z_score': z,
            'percentile': pct,
            'direction_match': dir_match,
            'cell_line': cell_line,
            'note': note,
            'per_cell_type': per_ct if len(per_ct) > 1 else None,
        }

    else:
        return {
            'grade': 'WEAK',
            'z_score': None,
            'percentile': None,
            'direction_match': None,
            'cell_line': None,
            'note': f'{upstream} was knocked down but {downstream} did not respond above threshold.',
        }


def validate_claims(claims_df, effects_df, stats_df=None):
    """Validate a batch of claims."""
    results = []
    for _, row in claims_df.iterrows():
        result = validate_single_claim(
            row['upstream_gene'],
            row['downstream_gene'],
            row['predicted_direction'],
            effects_df,
            stats_df,
        )
        result['upstream_gene'] = row['upstream_gene']
        result['downstream_gene'] = row['downstream_gene']
        result['predicted_direction'] = row['predicted_direction']
        results.append(result)
    return pd.DataFrame(results)


# ============================================================
# Streamlit App
# ============================================================

def main():
    st.set_page_config(
        page_title="TruthSeq",
        page_icon="🧬",
        layout="wide",
    )

    # Header
    st.title("TruthSeq")
    st.markdown(
        "Validate gene regulatory claims against real knockdown data from human cells. "
        "Upload your predictions and find out which ones hold up."
    )

    # Sidebar: data loading
    st.sidebar.header("Data")

    # Check for local data files
    data_files = find_data_files()

    effects_path = st.sidebar.text_input(
        "Effects parquet path",
        value=data_files.get('K562_effects', 'replogle_knockdown_effects.parquet'),
    )
    stats_path = st.sidebar.text_input(
        "Stats parquet path (optional)",
        value=data_files.get('K562_stats', 'replogle_knockdown_stats.parquet'),
    )

    # Check for RPE1
    rpe1_effects = data_files.get('RPE1_effects', '')
    if rpe1_effects:
        use_rpe1 = st.sidebar.checkbox("Also validate against RPE1 cells", value=True)
    else:
        use_rpe1 = False

    # Load data
    if not os.path.exists(effects_path):
        st.warning(
            "Knockdown data not found. Run `python3 setup.py` first to download "
            "the Replogle Perturb-seq data (~357 MB), then restart this app."
        )
        st.code("python3 setup.py\nstreamlit run app.py", language="bash")
        return

    effects_df, stats_df = load_replogle_data(effects_path, stats_path)

    # Optionally add RPE1
    if use_rpe1 and rpe1_effects:
        rpe1_eff, rpe1_stats = load_replogle_data(
            rpe1_effects,
            data_files.get('RPE1_stats', ''),
        )
        effects_df = pd.concat([effects_df, rpe1_eff], ignore_index=True)
        if rpe1_stats is not None and stats_df is not None:
            stats_df = pd.concat([stats_df, rpe1_stats], ignore_index=True)

    n_knockdowns = effects_df['knocked_down_gene'].nunique()
    n_pairs = len(effects_df)
    cell_types = effects_df['cell_line'].unique() if 'cell_line' in effects_df.columns else ['K562']

    st.sidebar.success(
        f"Loaded: {n_knockdowns:,} knockdowns, {n_pairs:,} gene pairs "
        f"({', '.join(cell_types)})"
    )

    # Main input area
    st.header("Your claims")

    input_method = st.radio(
        "How do you want to enter your predictions?",
        ["Upload CSV", "Enter manually", "Try example"],
        horizontal=True,
    )

    claims_df = None

    if input_method == "Upload CSV":
        uploaded = st.file_uploader(
            "Upload a CSV with columns: upstream_gene, downstream_gene, predicted_direction",
            type=['csv', 'tsv'],
        )
        if uploaded:
            sep = '\t' if uploaded.name.endswith('.tsv') else ','
            claims_df = pd.read_csv(uploaded, sep=sep)

    elif input_method == "Enter manually":
        st.markdown(
            "Enter one claim per line: `upstream_gene, downstream_gene, direction` "
            "(direction is UP or DOWN)"
        )
        text_input = st.text_area(
            "Claims",
            value="MYT1L, SCN2A, UP\nTCF4, CACNA1A, UP\nMEF2C, GRIN2B, UP",
            height=150,
        )
        if text_input.strip():
            lines = [l.strip() for l in text_input.strip().split('\n') if l.strip()]
            rows = []
            for line in lines:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 3:
                    rows.append({
                        'upstream_gene': parts[0],
                        'downstream_gene': parts[1],
                        'predicted_direction': parts[2].upper(),
                    })
            if rows:
                claims_df = pd.DataFrame(rows)

    elif input_method == "Try example":
        claims_df = pd.DataFrame([
            {'upstream_gene': 'SLC30A1', 'downstream_gene': 'MT2A', 'predicted_direction': 'DOWN'},
            {'upstream_gene': 'GATA1', 'downstream_gene': 'TYROBP', 'predicted_direction': 'DOWN'},
            {'upstream_gene': 'MYT1L', 'downstream_gene': 'SCN2A', 'predicted_direction': 'UP'},
            {'upstream_gene': 'TCF4', 'downstream_gene': 'CACNA1A', 'predicted_direction': 'UP'},
            {'upstream_gene': 'CHMP6', 'downstream_gene': 'SOD2', 'predicted_direction': 'DOWN'},
            {'upstream_gene': 'TP53', 'downstream_gene': 'CDKN1A', 'predicted_direction': 'UP'},
        ])
        st.info("Using example claims — includes known-true, known-weak, and untestable cases.")

    if claims_df is not None:
        # Validate required columns
        required = ['upstream_gene', 'downstream_gene', 'predicted_direction']
        missing = [c for c in required if c not in claims_df.columns]
        if missing:
            st.error(f"Missing required columns: {missing}")
            st.markdown("Your CSV needs these columns: `upstream_gene`, `downstream_gene`, `predicted_direction`")
            return

        claims_df['predicted_direction'] = claims_df['predicted_direction'].str.upper().str.strip()

        st.markdown(f"**{len(claims_df)} claims to validate**")

        # Run validation
        with st.spinner("Validating against knockdown data..."):
            results_df = validate_claims(claims_df, effects_df, stats_df)

        # Results summary
        st.header("Results")

        # Grade counts
        grade_counts = results_df['grade'].value_counts()
        cols = st.columns(5)
        for i, grade in enumerate(['VALIDATED', 'PARTIALLY_SUPPORTED', 'WEAK', 'CONTRADICTED', 'UNTESTABLE']):
            count = grade_counts.get(grade, 0)
            with cols[i]:
                st.metric(
                    label=grade.replace('_', ' ').title(),
                    value=count,
                )

        st.divider()

        # Per-claim results
        for _, row in results_df.iterrows():
            grade = row['grade']
            color = GRADE_COLORS.get(grade, '#95a5a6')

            col1, col2 = st.columns([3, 1])
            with col1:
                st.markdown(
                    f"**{row['upstream_gene']}** → **{row['downstream_gene']}** "
                    f"(predicted: {row['predicted_direction']})"
                )
                if row['note']:
                    st.caption(row['note'])
            with col2:
                st.markdown(
                    f"<div style='background-color:{color}; color:white; padding:8px 16px; "
                    f"border-radius:4px; text-align:center; font-weight:bold;'>"
                    f"{grade.replace('_', ' ')}</div>",
                    unsafe_allow_html=True,
                )

            # Per-cell-type details if available
            if row.get('per_cell_type'):
                with st.expander("Cell type breakdown"):
                    for ct_result in row['per_cell_type']:
                        st.markdown(
                            f"- **{ct_result['cell_line']}**: Z={ct_result['z_score']:.2f}, "
                            f"{ct_result['percentile']:.0f}th percentile, "
                            f"direction {'matches' if ct_result['direction_match'] else 'OPPOSES'}"
                        )

            st.divider()

        # Download results
        st.header("Download")
        csv_output = results_df.to_csv(index=False)
        st.download_button(
            "Download results as CSV",
            csv_output,
            file_name="truthseq_results.csv",
            mime="text/csv",
        )

        # Grade legend
        with st.expander("What do the grades mean?"):
            for grade, desc in GRADE_DESCRIPTIONS.items():
                color = GRADE_COLORS[grade]
                st.markdown(
                    f"<span style='color:{color}; font-weight:bold;'>{grade.replace('_', ' ')}</span>: {desc}",
                    unsafe_allow_html=True,
                )

    # Footer
    st.divider()
    st.caption(
        "TruthSeq validates gene regulatory predictions against the Replogle Perturb-seq atlas "
        "(Replogle et al. 2022, Cell). Data from Figshare. "
        "[GitHub](https://github.com/rsflinn/truthseq) · "
        "[About the method](https://github.com/rsflinn/truthseq#how-the-grades-work)"
    )


if __name__ == '__main__':
    main()
