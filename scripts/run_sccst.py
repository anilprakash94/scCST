#!/usr/bin/env python
"""Run scCST on an AnnData file.

Example:
    python scripts/run_sccst.py \
        --input data/processed/input.h5ad \
        --output data/processed/sccst_result.h5ad \
        --condition-col condition \
        --case Disease \
        --control Control \
        --batch-col batch \
        --cell-type-col Cell_Type \
        --latent-key X_pca_harmony
"""

import argparse
import scanpy as sc

from sccst import process_cell_type, make_results_adata


def parse_args():
    parser = argparse.ArgumentParser(description="Run scCST on an AnnData object.")
    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument("--output", required=True, help="Output .h5ad file")
    parser.add_argument("--condition-col", required=True)
    parser.add_argument("--case", required=True, help="Case/disease condition label")
    parser.add_argument("--control", required=True, help="Control condition label")
    parser.add_argument("--batch-col", required=True)
    parser.add_argument("--cell-type-col", required=True)
    parser.add_argument("--latent-key", default="X_pca_harmony")
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--n-neighborhoods", type=int, default=20)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--n-anchors-per-celltype", type=int, default=None)
    parser.add_argument("--top-genes-per-list", type=int, default=100)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)
    cell_importances, gene_names, run_info = process_cell_type(
        adata=adata,
        conditions=(args.case, args.control),
        condition_col=args.condition_col,
        batch_col=args.batch_col,
        sample_col=args.sample_col,
        n_neighborhoods=args.n_neighborhoods,
        cell_type_col=args.cell_type_col,
        n_jobs=args.n_jobs,
        latent_key=args.latent_key,
        top_genes_per_list=args.top_genes_per_list,
        fdr_threshold=args.fdr_threshold,
        n_anchors_per_celltype=args.n_anchors_per_celltype,
    )
    result = make_results_adata(adata, cell_importances, run_info, prefix="rra", subset_to_results=True)
    result.write(args.output)


if __name__ == "__main__":
    main()
