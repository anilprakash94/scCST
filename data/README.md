# Data

Large datasets and generated result files should not be committed to GitHub.

Recommended local layout:

```text
data/
├── raw/
├── processed/
└── external/
```

Use this directory only for instructions, metadata, and small toy files. For manuscript-sized `.h5ad`, `.h5`, `.loom`, `.rds`, or result files, use an external archive such as Zenodo, Figshare, OSF, Dryad, institutional storage, or a GitHub Release if the files are small enough.
