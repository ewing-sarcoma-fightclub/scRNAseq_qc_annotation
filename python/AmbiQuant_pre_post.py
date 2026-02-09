#!/usr/bin/env python3
import argparse
import csv
import os
import sys
from pathlib import Path

import scanpy as sc
from scipy.io import mmread

SCRIPT_DIR = Path(__file__).resolve().parent
abq = None


def _load_ambiquant(ambiquant_root):
    global abq
    amb_root = Path(ambiquant_root)
    if not amb_root.exists():
        raise SystemExit(f"AmbiQuant repo not found at {amb_root}")
    sys.path.insert(0, str(amb_root))
    import AmbiQuantFunctions as abq_local  # noqa: E402

    abq = abq_local


def _write_metrics(out_csv, sample, stage, ret, overall_score):
    fields = [
        "sample",
        "stage",
        "empty_droplet_area_ratio",
        "inv_max_secant",
        "inv_secant_std",
        "inv_auc_cum_curve",
        "num_amb_genes",
        "mean_pct_ambient",
        "overall_score",
    ]
    row = {
        "sample": sample,
        "stage": stage,
        "empty_droplet_area_ratio": ret[0],
        "inv_max_secant": ret[1],
        "inv_secant_std": ret[2],
        "inv_auc_cum_curve": ret[3],
        "num_amb_genes": ret[4],
        "mean_pct_ambient": ret[5],
        "overall_score": overall_score,
    }

    file_exists = os.path.exists(out_csv)
    with open(out_csv, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def run_pre(raw_dir, out_dir, sample, mito_tag, inflection_fold, max_cell, dropout_thresh):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_10x_mtx(str(raw_dir), var_names="gene_symbols", cache=False)
    dat = abq.cut_off_h5ad(
        adata,
        inflection_fold=inflection_fold,
        est_real_cell=None,
        max_cell=max_cell,
        qc_pipe_inflection=True,
        mito_tag=mito_tag,
        run_qc=True,
        print_info=True,
    )

    fig_path = out_dir / f"{sample}_pre_AmbiQuant.png"
    amb_path = out_dir / f"{sample}_pre_ambient_genes.txt"

    ret = abq.formatted_figures_inverted(
        dat,
        save_amb_ls=str(amb_path),
        save_fig=str(fig_path),
        show_dat_name=f"{sample} (pre-SoupX)",
        invert_scores=True,
        dropout_thresh=dropout_thresh,
        ncols=2,
    )

    overall = abq.calc.overall_score(ret, from_formatted_figures=True)
    _write_metrics(out_dir / "AmbiQuant_metrics.csv", sample, "pre", ret, overall)


def run_post(filtered_dir, soupx_mtx, out_dir, sample, mito_tag, inflection_fold, max_cell, dropout_thresh):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_10x_mtx(str(filtered_dir), var_names="gene_symbols", cache=False)
    soupx = mmread(str(soupx_mtx)).tocsr()
    if soupx.shape != adata.X.shape:
        soupx = soupx.T.tocsr()
    if soupx.shape != adata.X.shape:
        raise ValueError(
            f"SoupX matrix shape {soupx.shape} does not match filtered matrix {adata.X.shape}"
        )

    adata.X = soupx

    dat = abq.cut_off_h5ad(
        adata,
        inflection_fold=inflection_fold,
        est_real_cell=None,
        max_cell=max_cell,
        qc_pipe_inflection=True,
        mito_tag=mito_tag,
        run_qc=True,
        print_info=True,
    )

    fig_path = out_dir / f"{sample}_post_AmbiQuant.png"
    amb_path = out_dir / f"{sample}_post_ambient_genes.txt"

    ret = abq.formatted_figures_inverted(
        dat,
        save_amb_ls=str(amb_path),
        save_fig=str(fig_path),
        show_dat_name=f"{sample} (post-SoupX)",
        invert_scores=True,
        dropout_thresh=dropout_thresh,
        ncols=2,
    )

    overall = abq.calc.overall_score(ret, from_formatted_figures=True)
    _write_metrics(out_dir / "AmbiQuant_metrics.csv", sample, "post", ret, overall)


def main():
    p = argparse.ArgumentParser(description="Run AmbiQuant before and after SoupX.")
    p.add_argument("--root", required=True, help="Root directory with per-sample subfolders")
    p.add_argument("--out", required=True, help="Output directory for AmbiQuant results")
    p.add_argument("--mito-tag", default="MT-", help="Mitochondrial gene prefix (default: MT-)")
    p.add_argument("--inflection-fold", type=int, default=4)
    p.add_argument("--max-cell", type=int, default=10000)
    p.add_argument("--dropout-thresh", type=float, default=2)
    p.add_argument("--soupx-suffix", default="_soupx.mtx")
    p.add_argument("--soupx-root", default=None)
    p.add_argument("--ambiquant-root", default=os.environ.get("AMBIQUANT_REPO"))
    p.add_argument("--skip-pre", action="store_true")
    p.add_argument("--skip-post", action="store_true")
    args = p.parse_args()

    if not args.ambiquant_root:
        raise SystemExit("AmbiQuant repo path not provided. Set --ambiquant-root or AMBIQUANT_REPO.")
    _load_ambiquant(args.ambiquant_root)

    root = Path(args.root)
    out = Path(args.out)

    for sample_dir in sorted(root.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name

        raw_dir = sample_dir / "raw_feature_bc_matrix"
        filtered_dir = sample_dir / "filtered_feature_bc_matrix"
        if not raw_dir.exists():
            raw_dir = sample_dir / "outs" / "raw_feature_bc_matrix"
        if not filtered_dir.exists():
            filtered_dir = sample_dir / "outs" / "filtered_feature_bc_matrix"
        if args.soupx_root:
            soupx_mtx = Path(args.soupx_root) / sample / f"{sample}{args.soupx_suffix}"
        else:
            soupx_mtx = sample_dir / "SoupX" / f"{sample}{args.soupx_suffix}"

        if not args.skip_pre:
            if raw_dir.exists():
                run_pre(raw_dir, out, sample, args.mito_tag, args.inflection_fold, args.max_cell, args.dropout_thresh)
            else:
                print(f"[WARN] Missing raw_feature_bc_matrix for {sample}; skipping pre-SoupX AmbiQuant")

        if not args.skip_post:
            if filtered_dir.exists() and soupx_mtx.exists():
                run_post(filtered_dir, soupx_mtx, out, sample, args.mito_tag, args.inflection_fold, args.max_cell, args.dropout_thresh)
            else:
                print(f"[WARN] Missing filtered_feature_bc_matrix or SoupX matrix for {sample}; skipping post-SoupX AmbiQuant")


if __name__ == "__main__":
    main()
