# Single-core runtime benchmark

Alternative versions of all supplied benchmark scripts.

- Five datasets only: base cell counts 2,000, 5,000, 10,000, 15,000, and 20,000.
- No replicates.
- Fixed defaults: responder percent 40, dropout 0.4, Harmony theta 2.
- One core at every level.
- Common neighborhood size 20.

Benchmark totals:
- scCST = preprocessing + `process_cell_type`.
- miloDE = preprocessing + `de_test_neighbourhoods`. Assignment and AUC setup are reported separately, not added to the benchmark total.
- CACOA = preprocessing + `estimateClusterFreeDE`. Neighbor finding/object construction are reported separately, not added to the benchmark total.
- Output writing is excluded.

Run in this order: simulation, conversion, both Stage-1 preprocessing scripts, both Stage-2 DE scripts, scCST, then the summary combiner.
