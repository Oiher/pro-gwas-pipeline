#!/usr/bin/env python3
"""
Manhattan and QQ plot generator for GWAS results

Usage:
    manhattan.py --input <results.tsv> --model <glm|cph|lmm_gallop> [--suffix <str>]

Arguments:
    --input: Path to GWAS results TSV file
    --model: Analysis model type (glm, cph, or lmm_gallop)
    --suffix: Optional suffix for cohort naming (default: "null")
"""

import argparse
import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
from qmplot import manhattanplot, qqplot

# Significance thresholds for Manhattan plot labeling.
# NOMINAL_P: "suggestive"/nominal significance -- also used as qmplot's
#   sign_marker_p, which is what actually determines which variants become
#   label candidates (qmplot's suggestiveline/genomewideline only draw the
#   horizontal reference lines, they don't drive labeling). Set to the nominal
#   threshold rather than MTC_P so BOTH nominally-significant and
#   MTC-significant variants get labeled -- MTC significance is a strict
#   subset of nominal significance, so this single threshold covers both.
# MTC_P: standard GWAS genome-wide (multiple-testing-corrected) significance
#   convention, used only for the reference line -- not recomputed per-run
#   from the actual number of variants tested.
NOMINAL_P = 1e-5
MTC_P = 5e-8


def plot_summary_stats(data, cohort, outcome, model):
    """Generate Manhattan and QQ plots for GWAS summary statistics.

    Args:
        data: pandas DataFrame with GWAS results
        cohort: Cohort name for plot titles and filenames
        outcome: Phenotype/outcome name
        model: Model type (glm, cph, or lmm_gallop)
    """
    xtick = set(['chr' + i for i in list(map(str, range(1, 14))) + ['15', '17', '19', '22']])

    # Shared kwargs: draws both significance lines, colors/labels variants
    # crossing NOMINAL_P (which includes anything crossing MTC_P too) with
    # their variant ID -- qmplot dedupes to one label per ~50kb LD block
    # (ld_block_size) so a cluster of correlated significant SNPs at the same
    # locus doesn't produce overlapping labels.
    sig_kws = dict(
        snp="ID",
        suggestiveline=NOMINAL_P,
        genomewideline=MTC_P,
        sign_marker_p=NOMINAL_P,
        is_annotate_topsnp=True,
    )

    if model == "lmm_gallop":
        # GALLOP longitudinal model - plot both intercept and slope

        # Intercept Manhattan plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan Intercept {cohort} {outcome}",
                      pv="Pi", ax=ax,
                      xtick_label_set=xtick,
                      **sig_kws)
        plt.savefig(f"{cohort}_{outcome}_manhattan_intercept.{model}.png", dpi=300)

        # Intercept QQ plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        qqplot(data=data["Pi"],
               marker="o",
               title=f"QQ Intercept {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq_intercept.{model}.png", dpi=300)

        # Slope Manhattan plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan Slope {cohort} {outcome}",
                      pv="Ps", ax=ax,
                      xtick_label_set=xtick,
                      **sig_kws)
        plt.savefig(f"{cohort}_{outcome}_manhattan_slope.{model}.png", dpi=300)
        
        # Slope QQ plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        qqplot(data=data["Ps"],
               marker="o",
               title=f"QQ Slope {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq_slope.{model}.png", dpi=300)
    
    elif model in ["glm", "cph"]:
        # Standard GLM or Cox PH model - single P-value
        
        # Manhattan plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan {model} {cohort} {outcome}",
                      pv="P", ax=ax,
                      xtick_label_set=xtick,
                      **sig_kws)
        plt.savefig(f"{cohort}_{outcome}_manhattan.{model}.png", dpi=300)
        
        # QQ plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        qqplot(data=data["P"],
               marker="o",
               title=f"QQ {model} {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq.{model}.png", dpi=300)
    
    else:
        raise ValueError(f"Model '{model}' not recognized. Use 'glm', 'cph', or 'lmm_gallop'.")
    
    print(f"✅ Results plotting success for {cohort}_{outcome} ({model})")


def main():
    """Main entry point for Manhattan plot generation."""
    parser = argparse.ArgumentParser(
        description="Generate Manhattan and QQ plots from GWAS results",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", required=True, help="Input TSV file with GWAS results")
    parser.add_argument("--model", required=True, choices=["glm", "cph", "lmm_gallop"],
                       help="Analysis model type")
    parser.add_argument("--suffix", default="null", help="Cohort suffix (default: null)")
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"❌ Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Extract metadata from filename
    # Expected format: {ancestry}_{cohort}_{phenotype}_allresults.tsv
    gwas_name = os.path.splitext(os.path.basename(args.input))[0]
    
    try:
        cohort = gwas_name.split('_')[1]
        find_phenoname = cohort + r'_(.*?)_allresults'
        result = re.search(find_phenoname, gwas_name)
        if not result:
            raise ValueError(f"Could not extract phenotype name from filename: {gwas_name}")
        pheno = result.group(1)
    except (IndexError, ValueError) as e:
        print(f"❌ Error parsing filename: {e}", file=sys.stderr)
        print(f"Expected format: <ancestry>_<cohort>_<phenotype>_allresults.tsv", file=sys.stderr)
        sys.exit(1)
    
    # Add suffix to cohort if provided
    cohort_suffix = f"{cohort}.{args.suffix}" if args.suffix and args.suffix != "null" else cohort
    
    # Load and prepare data
    print(f"📊 Loading data from: {args.input}")
    try:
        df = pd.read_csv(args.input, sep="\t")
    except Exception as e:
        print(f"❌ Error reading file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"   Loaded {len(df)} variants")
    
    # Clean data
    df = df.dropna(how="any", axis=0)
    print(f"   After removing missing values: {len(df)} variants")
    
    # Add chromosome ordering for proper sorting
    df['chr_order'] = df['#CHROM'].str.replace('chr', '')
    df['chr_order'] = df['chr_order'].astype(int)
    df = df.sort_values(by=['chr_order', 'POS'])
    
    # Generate plots
    print(f"🎨 Generating plots for cohort={cohort_suffix}, phenotype={pheno}, model={args.model}")
    plot_summary_stats(data=df, cohort=cohort_suffix, outcome=pheno, model=args.model)


if __name__ == "__main__":
    main()
