#!/usr/bin/env python
"""Minimal downstream scCST analysis helper."""

import argparse
import scanpy as sc

from sccst.downstream import cluster_on_gene_scores, analyze_global_drivers


def parse_args():
    parser = argparse.ArgumentParser(description="Run lightweight downstream analysis on scCST output.")
    parser.add_argument("--input", required=True, help="Input scCST result .h5ad file")
    parser.add_argument("--output", required=True, help="Output .h5ad file")
    parser.add_argument("--n-top-genes", type=int, default=2000)
    return parser.parse_args()


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)
    adata = cluster_on_gene_scores(adata, n_top_genes=args.n_top_genes)
    drivers = analyze_global_drivers(adata)
    drivers.to_csv(args.output.replace(".h5ad", "_global_drivers.csv"), index=False)
    adata.write(args.output)


if __name__ == "__main__":
    main()
