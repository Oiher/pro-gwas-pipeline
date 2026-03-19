#!/usr/bin/env python3
"""
Forest plot generator for stratified survival GWAS results.

For each unique variant found across all files with the given extension, produces
a forest plot showing per-study beta (log HR) and 95% CI, allele frequency, and
sample size. Studies whose |SE| exceeds SE_THRESHOLD are shown as a diamond but
are excluded from the x-axis scale (so runaway CIs from rare-variant convergence
failures don't squash the readable studies).

Output: <results_dir>/forest_plot/<CHROM>_<POS>_<REF>_<ALT>.png

Usage:
    python make_forest_plot.py [results_dir] [extension]

Defaults:
    results_dir : proj_r11/analyses/focus_strata/surv/results
    extension   : coxph
"""

import os
import sys
import glob
import re

import math

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sp_stats

# ── Parameters ────────────────────────────────────────────────────────────────
RESULTS_DIR  = sys.argv[1] if len(sys.argv) > 1 else "proj_r11/analyses/focus_strata/surv/results"
EXTENSION    = sys.argv[2].lstrip(".") if len(sys.argv) > 2 else "coxph"
OUTPUT_DIR   = os.path.join(RESULTS_DIR, "forest_plot")
SE_THRESHOLD = 100        # |SE| > this → "extreme": point plotted but CI excluded from scale
DPI          = 150
COL_WIDTHS   = [3.5, 4.5, 3.8]   # [labels, forest, stats] panel widths (inches)
ROW_H        = 0.48               # figure height per study row (inches)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load all studies ──────────────────────────────────────────────────────────
print(f"Reading .{EXTENSION} files from {RESULTS_DIR} …", flush=True)
frames = []
for fpath in sorted(glob.glob(os.path.join(RESULTS_DIR, f"*.{EXTENSION}"))):
    study = os.path.basename(fpath).split("_EUR_")[0]
    try:
        df = pd.read_csv(fpath, sep="\t")
        df.columns = df.columns.str.lstrip("#").str.strip()
        df["STUDY"] = study
        frames.append(df)
    except Exception as exc:
        print(f"  [warn] {os.path.basename(fpath)}: {exc}", flush=True)

if not frames:
    sys.exit(f"No .{EXTENSION} files found in " + RESULTS_DIR)

df_all = pd.concat(frames, ignore_index=True)
for col in ("BETA", "SE", "A1_FREQ", "OBS_CT", "P"):
    df_all[col] = pd.to_numeric(df_all[col], errors="coerce")

print(f"  Loaded {len(frames)} studies, {len(df_all)} total rows", flush=True)


def _safe(s: str) -> str:
    """Replace characters unsafe for filenames with underscores."""
    return re.sub(r"[^\w.-]", "_", str(s))


def compute_meta(betas, ses):
    """
    Inverse-variance weighted fixed-effect and DerSimonian-Laird random-effect
    meta-analysis.

    Parameters
    ----------
    betas, ses : array-like (non-extreme studies only)

    Returns
    -------
    dict with fe_*, re_*, tau2, Q, Q_p, k
    """
    betas = np.asarray(betas, dtype=float)
    ses   = np.asarray(ses,   dtype=float)
    # drop studies with zero or non-finite SE before computing weights
    valid = (ses > 0) & np.isfinite(ses) & np.isfinite(betas)
    betas = betas[valid]
    ses   = ses[valid]
    k     = len(betas)
    if k == 0:
        return None

    w = 1.0 / ses**2
    W = w.sum()

    # Fixed effect (IVW)
    fe_beta = (w * betas).sum() / W
    fe_se   = math.sqrt(1.0 / W)
    fe_z    = fe_beta / fe_se
    fe_p    = 2.0 * sp_stats.norm.sf(abs(fe_z))

    # Cochran's Q
    Q   = float((w * (betas - fe_beta)**2).sum())
    Q_p = float(sp_stats.chi2.sf(Q, df=k - 1)) if k > 1 else float("nan")

    # DerSimonian-Laird tau^2 (undefined / zero for k=1)
    if k > 1:
        C    = W - (w**2).sum() / W
        tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    else:
        tau2 = 0.0

    # Random effect (equals fixed effect when tau2=0 or k=1)
    w_re    = 1.0 / (ses**2 + tau2)
    re_beta = (w_re * betas).sum() / w_re.sum()
    re_se   = math.sqrt(1.0 / w_re.sum())
    re_z    = re_beta / re_se
    re_p    = 2.0 * sp_stats.norm.sf(abs(re_z))

    # I² heterogeneity statistic
    i2 = max(0.0, (Q - (k - 1)) / Q * 100.0) if (k > 1 and Q > 0) else 0.0

    return dict(
        fe_beta=fe_beta, fe_se=fe_se,
        fe_ci_lo=fe_beta - 1.96*fe_se, fe_ci_hi=fe_beta + 1.96*fe_se, fe_p=fe_p,
        re_beta=re_beta, re_se=re_se,
        re_ci_lo=re_beta - 1.96*re_se, re_ci_hi=re_beta + 1.96*re_se, re_p=re_p,
        tau2=tau2, Q=Q, Q_p=Q_p, i2=i2, k=k,
    )


def _draw_meta_row(ax_l, ax_m, ax_r, y, label, prefix, meta, x_lo, x_hi, color):
    """Draw one meta-analysis row (fixed or random effect)."""
    beta  = meta[f"{prefix}_beta"]
    ci_lo = meta[f"{prefix}_ci_lo"]
    ci_hi = meta[f"{prefix}_ci_hi"]
    p     = meta[f"{prefix}_p"]

    ax_l.text(0.02, y, label, va="center", fontsize=8.5, fontweight="bold", color=color)

    CAP = 0.13
    ax_m.plot([max(ci_lo, x_lo), min(ci_hi, x_hi)], [y, y],
              color=color, linewidth=2.5, zorder=2, solid_capstyle="butt")
    if ci_lo >= x_lo:
        ax_m.plot([ci_lo, ci_lo], [y - CAP, y + CAP], color=color, linewidth=2.5)
    if ci_hi <= x_hi:
        ax_m.plot([ci_hi, ci_hi], [y - CAP, y + CAP], color=color, linewidth=2.5)
    ax_m.scatter([beta], [y], marker="D", s=72, color=color, zorder=4, edgecolors="none")

    b_str = f"{beta:.3f} ({ci_lo:.3f}, {ci_hi:.3f})"
    p_str = (f"{p:.2e}" if p < 0.001 else f"{p:.3f}") if not math.isnan(p) else "NA"
    ax_r.text(0.02, y, b_str,  va="center", fontsize=7,   color=color)
    ax_r.text(0.98, y, p_str,  va="center", ha="right", fontsize=7.5, color=color)


# ── Meta row y positions (below cohort rows, separated at META_SEP_Y) ─────────
META_SEP_Y = 0.4
Y_FE       = -0.4
Y_RE       = -1.4
Y_LO_BASE  = -2.2

# ── Output TSV collectors ─────────────────────────────────────────────────────
study_rows = []   # one row per (variant × study)
meta_rows  = []   # one row per variant (meta-analysis summary)

# ── One forest plot per variant ───────────────────────────────────────────────
n_saved = 0

for (chrom, pos, ref, alt), grp in df_all.groupby(["CHROM", "POS", "REF", "ALT"], sort=False):
    grp = (grp.dropna(subset=["BETA", "SE"])
              .sort_values("STUDY")
              .reset_index(drop=True))
    if grp.empty:
        continue

    n = len(grp)
    grp["extreme"] = grp["SE"].abs() > SE_THRESHOLD
    grp["ci_lo"]   = grp["BETA"] - 1.96 * grp["SE"]
    grp["ci_hi"]   = grp["BETA"] + 1.96 * grp["SE"]

    # Non-extreme studies
    normal = grp[~grp["extreme"]]
    if normal.empty:
        continue  # skip variant if every study has extreme SE

    # Meta-analysis (fixed-effect + DerSimonian-Laird random-effect)
    meta = compute_meta(normal["BETA"].values, normal["SE"].values)
    if meta is None:
        continue

    # ── Accumulate per-study rows ─────────────────────────────────────────────
    for _, row in grp.iterrows():
        study_rows.append(dict(
            CHROM=chrom, POS=pos, REF=ref, ALT=alt,
            STUDY=row["STUDY"],
            BETA=row["BETA"], SE=row["SE"],
            CI_LO=row["BETA"] - 1.96 * row["SE"],
            CI_HI=row["BETA"] + 1.96 * row["SE"],
            P=row["P"],
            A1_FREQ=row["A1_FREQ"],
            OBS_CT=row["OBS_CT"],
            EXTREME=bool(row["extreme"]),
        ))

    # ── Accumulate meta-analysis summary row ──────────────────────────────────
    af_vals   = normal["A1_FREQ"].dropna()
    total_n   = normal["OBS_CT"].dropna().sum()
    meta_rows.append(dict(
        CHROM=chrom, POS=pos, REF=ref, ALT=alt,
        N_STUDIES=meta["k"],
        N_STUDIES_EXTREME=int(grp["extreme"].sum()),
        TOTAL_N=int(total_n),
        A1_FREQ_MEAN=af_vals.mean() if len(af_vals) else float("nan"),
        A1_FREQ_SD=af_vals.std()    if len(af_vals) > 1 else float("nan"),
        FE_BETA=meta["fe_beta"],   FE_SE=meta["fe_se"],
        FE_CI_LO=meta["fe_ci_lo"], FE_CI_HI=meta["fe_ci_hi"], FE_P=meta["fe_p"],
        RE_BETA=meta["re_beta"],   RE_SE=meta["re_se"],
        RE_CI_LO=meta["re_ci_lo"], RE_CI_HI=meta["re_ci_hi"], RE_P=meta["re_p"],
        I2=meta["i2"], Q=meta["Q"], Q_P=meta["Q_p"], TAU2=meta["tau2"],
    ))

    # x-axis range: non-extreme study CIs + meta CIs
    x_vals = (list(normal["ci_lo"]) + list(normal["ci_hi"]) +
              [meta["fe_ci_lo"], meta["fe_ci_hi"], meta["re_ci_lo"], meta["re_ci_hi"]])
    x_lo = min(x_vals)
    x_hi = max(x_vals)
    if np.isclose(x_lo, x_hi):
        x_lo -= 1.0
        x_hi += 1.0
    pad  = (x_hi - x_lo) * 0.12
    x_lo -= pad
    x_hi += pad

    # y layout: cohort rows n→1, meta rows below separator
    ys   = np.arange(n, 0, -1, dtype=float)
    y_lo = Y_LO_BASE
    y_hi = n + 0.75

    fig_h = max(4.5, (n + 3) * ROW_H + 1.8)
    fig, (ax_l, ax_m, ax_r) = plt.subplots(
        1, 3,
        figsize=(sum(COL_WIDTHS), fig_h),
        gridspec_kw={"width_ratios": COL_WIDTHS},
    )

    # ── Setup axes ────────────────────────────────────────────────────────────
    for ax in (ax_l, ax_r):
        ax.set_xlim(0, 1)
        ax.set_ylim(y_lo, y_hi)
        ax.axis("off")

    ax_m.set_xlim(x_lo, x_hi)
    ax_m.set_ylim(y_lo, y_hi)
    ax_m.spines[["top", "right", "left"]].set_visible(False)
    ax_m.set_yticks([])
    ax_m.tick_params(left=False)
    ax_m.set_xlabel("Beta (log HR)", fontsize=8.5)
    ax_m.set_title(f"{chrom}:{pos}  {ref} → {alt}", fontsize=9.5, fontweight="bold", pad=6)
    ax_m.axvline(0, color="#555", linestyle="--", linewidth=0.8, zorder=1)

    # ── Header row ────────────────────────────────────────────────────────────
    hdr_y = n + 0.56
    for ax in (ax_l, ax_m, ax_r):
        ax.axhline(n + 0.28, color="#bbb", linewidth=0.6)

    ax_l.text(0.02, hdr_y, "Study",      fontweight="bold", fontsize=8.5, va="center")
    ax_l.text(0.98, hdr_y, "AF / N",     fontweight="bold", fontsize=8.5, va="center", ha="right")
    ax_r.text(0.02, hdr_y, "β (95% CI)", fontweight="bold", fontsize=8.5, va="center")
    ax_r.text(0.98, hdr_y, "P",          fontweight="bold", fontsize=8.5, va="center", ha="right")

    # ── Data rows ─────────────────────────────────────────────────────────────
    for i, row in grp.iterrows():
        y         = ys[i]
        ext       = bool(row["extreme"])
        pt_color  = "steelblue" if not ext else "lightgray"
        txt_color = "black"    if not ext else "dimgray"
        beta      = float(row["BETA"])

        # Left panel: study name  |  AF / N
        ax_l.text(0.02, y, row["STUDY"], va="center", fontsize=8, color=txt_color)
        af_str = f"{row['A1_FREQ']:.3f}" if pd.notna(row["A1_FREQ"]) else "NA"
        n_str  = str(int(row["OBS_CT"])) if pd.notna(row["OBS_CT"])  else "NA"
        ax_l.text(0.98, y, f"{af_str} / {n_str}",
                  va="center", ha="right", fontsize=7.5, color="#666")

        # Middle panel: CI line + point estimate
        if ext:
            # show only the point estimate as a grey diamond; no CI bar
            ax_m.scatter([np.clip(beta, x_lo, x_hi)], [y],
                         marker="D", s=36, color=pt_color, alpha=0.7, zorder=3)
        else:
            lo, hi = float(row["ci_lo"]), float(row["ci_hi"])
            CAP = 0.10
            # CI line clipped to plot range
            ax_m.plot([max(lo, x_lo), min(hi, x_hi)], [y, y],
                      color=pt_color, linewidth=1.8, zorder=2, solid_capstyle="butt")
            # end caps (only drawn when the CI end is within range)
            if lo >= x_lo:
                ax_m.plot([lo, lo], [y - CAP, y + CAP], color=pt_color, linewidth=1.8)
            if hi <= x_hi:
                ax_m.plot([hi, hi], [y - CAP, y + CAP], color=pt_color, linewidth=1.8)
            # point estimate square
            ax_m.scatter([beta], [y], marker="s", s=52,
                         color=pt_color, zorder=4, edgecolors="none")

        # Right panel: numeric summary
        p = row["P"]
        if ext:
            b_str = f"{beta:.2f}  [extreme SE]"
        else:
            b_str = f"{beta:.3f} ({float(row['ci_lo']):.3f}, {float(row['ci_hi']):.3f})"
        p_str = ("NA"        if pd.isna(p)   else
                 f"{p:.2e}"  if p < 0.001    else
                 f"{p:.3f}")
        ax_r.text(0.02, y, b_str,  va="center", fontsize=7,   color=txt_color)
        ax_r.text(0.98, y, p_str,  va="center", ha="right", fontsize=7.5, color=txt_color)

    # ── Meta-analysis separator ───────────────────────────────────────────────
    for ax in (ax_l, ax_m, ax_r):
        ax.axhline(META_SEP_Y, color="#888", linewidth=1.2, linestyle="-")

    q_str   = (f"Q={meta['Q']:.2f} (p={meta['Q_p']:.3f})"
               if not math.isnan(meta["Q_p"]) else "")
    tau_str = f"\u03c4\u00b2={meta['tau2']:.4f}"
    i2_str  = f"I\u00b2={meta['i2']:.1f}%"
    ax_l.text(0.02, META_SEP_Y + 0.23,
              f"Meta-analysis  (k={meta['k']}   {q_str}   {i2_str}   {tau_str})",
              va="center", fontsize=7.5, fontstyle="italic", color="#555")

    # ── Meta rows ─────────────────────────────────────────────────────────────
    _draw_meta_row(ax_l, ax_m, ax_r, Y_FE, "Fixed effect",  "fe", meta, x_lo, x_hi, "#c0392b")
    _draw_meta_row(ax_l, ax_m, ax_r, Y_RE, "Random effect", "re", meta, x_lo, x_hi, "#27ae60")

    # ── Footnote for extreme studies ──────────────────────────────────────────
    n_ext = int(grp["extreme"].sum())
    if n_ext:
        fig.text(
            0.5, 0.005,
            f"* {n_ext} {'study' if n_ext == 1 else 'studies'} with |SE| > {SE_THRESHOLD}"
            " shown as ◆ but excluded from x-axis scale and meta-analysis",
            ha="center", fontsize=7, color="gray", style="italic",
        )

    plt.tight_layout(pad=0.5)

    fname = f"{_safe(chrom)}_{pos}_{_safe(ref)}_{_safe(alt)}.png"
    plt.savefig(os.path.join(OUTPUT_DIR, fname), dpi=DPI, bbox_inches="tight")
    plt.close()
    n_saved += 1
    if n_saved % 100 == 0:
        print(f"  {n_saved} variants …", flush=True)

print(f"\nDone — {n_saved} forest plots saved to {OUTPUT_DIR}/")

# ── Write TSV outputs ─────────────────────────────────────────────────────────
study_tsv = os.path.join(OUTPUT_DIR, "per_study_results.tsv")
meta_tsv  = os.path.join(OUTPUT_DIR, "meta_results.tsv")

pd.DataFrame(study_rows).to_csv(study_tsv, sep="\t", index=False,
    float_format="%.6g")
print(f"Per-study results  → {study_tsv}")

df_meta = pd.DataFrame(meta_rows)
# Round display columns for readability
for col in ("A1_FREQ_MEAN", "A1_FREQ_SD",
            "FE_BETA", "FE_SE", "FE_CI_LO", "FE_CI_HI",
            "RE_BETA", "RE_SE", "RE_CI_LO", "RE_CI_HI",
            "I2", "Q", "TAU2"):
    df_meta[col] = df_meta[col].round(6)
for col in ("FE_P", "RE_P", "Q_P"):
    df_meta[col] = df_meta[col].apply(
        lambda x: f"{x:.4e}" if pd.notna(x) and float(x) < 0.001 else round(float(x), 6) if pd.notna(x) else "NA"
    )
df_meta.to_csv(meta_tsv, sep="\t", index=False)
print(f"Meta-analysis results → {meta_tsv}")
