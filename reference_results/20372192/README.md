# Reference experiment 20372192

This directory contains the numerical outputs of the completed CSU HPC run.
There are 24 datasets, 360 algorithm/partition groups, 360 measured calls,
and 360 separate warm-up calls. All recorded calls have status `OK`.
Every group has one measured repetition; times are not averages over repeated runs.

| File | Meaning |
|---|---|
| `plan.csv` | Original expanded experiment plan. |
| `results.tsv` | Byte-for-byte copy of `results-partial.tsv`, with 17-digit Kemeny/reference values. |
| `summary.csv` | Original per-group final summary. |
| `warmup-results.tsv` | Original warm-up observations. |
| `timing-matlab.tsv` | Original stage timings; TOTAL rows include nested child stages. |
| `provenance.json` | SHA256 hashes tying the release to the original source/data snapshot. |

The TSV preserves more printed Kemeny digits than the original final CSV export.
Use the recorded `RelativeError` column or full-precision reference values rather
than recomputing tiny differences from a table rounded for publication.

## Environment

- MATLAB 9.11.0.1769968 (R2021b), x86-64 Linux.
- Intel Xeon Gold 6248R, two sockets with 24 physical cores per socket;
  nominal frequency 3.00 GHz.
- Job allocation: 15 CPU cores and 100 GiB RAM, not the whole node's RAM capacity.
- Process pool: 14 workers, one computational thread per worker, reused throughout.
- Direct/Bini client thread limit: 15. Du/New client thread limit: 1.
- BLAS: Intel MKL 2019.0.3; LAPACK 3.7.0 reported by MATLAB.
- Full shell elapsed time: 69807.343372 seconds (19:23:27).
- The effective scheduler time limit is not known. An 8-hour directive in the
  separately downloaded submission script conflicts with the completed duration.

## Interpretation

The runner uses the mean of successful measured Direct results in each test row
as the numerical reference. Because there is only one measured Direct call here,
that mean equals its single result. This reference is not an exact solution.

For tables choosing the `k` that minimizes New's time, use the **same chosen k**
for Du's time and error. Do not combine Du's independently best k into that row.
Single-repetition minima are exploratory selections and need repeated measurements
before claiming stable performance rankings.

Only numeric result files are distributed. Raw machine/account paths, MATLAB
preferences, private certificates, and worker caches are excluded. The input
matrices and baseline licenses are documented in `../../matrix/README.md` and
`../../THIRD_PARTY_NOTICES.md`.
