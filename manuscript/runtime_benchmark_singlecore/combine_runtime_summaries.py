#!/usr/bin/env python
from pathlib import Path
import pandas as pd, os
base=Path(os.environ.get("RUNTIME_DATA_DIR", "/home/anilprakash/labs/Mei/projects/anil/srda/notebooks/data/scrna_seq/simulation"))
patterns=[("scCST",base/"sccst_runtime_5_cell_counts"/"sccst_runtime_summary.csv"),("miloDE",base/"milode_runtime_5_cell_counts_results"/"milode_hvg_rds_run_summary.csv"),("CACOA",base/"cacoa_runtime_5_cell_counts_results"/"cacoa_parameter_sweeps_run_summary_from_hvg_rds.csv")]
parts=[]
for method,p in patterns:
    if p.exists():
        d=pd.read_csv(p); d["method"]=method; parts.append(d)
    else: print("Missing",p)
out=pd.concat(parts,ignore_index=True,sort=False) if parts else pd.DataFrame()
out.to_csv(base/"runtime_benchmark_three_methods_singlecore.csv",index=False)
print(out)
