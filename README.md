# Kemeny-and-Gamma

MATLAB implementation of a gamma/feedback method for computing Kemeny's constant,
with reproducible comparisons against Direct, Bini-Rec, Bini-DirRec, and Du.
The proposed implementation is named **New-full-rhs** in the benchmark output.

The default configuration reproduces the experiment plan of HPC job `20372192`:
24 datasets, six partition counts (`4, 6, 8, 10, 12, 14`) for Du and New,
one full warm-up and one measured call for each enabled algorithm/partition group.
There are **360 measured calls and 360 warm-up calls**.
This full batch took about **19 hours 23 minutes** on the reference system.
Start with the small example below before running the complete batch.

## Requirements

- MATLAB with `dissect` and Parallel Computing Toolbox for the parallel methods.
- Reference experiments used MATLAB **9.11.0.1769968 (R2021b)**.
- The release checks use independent small matrices and were run with R2025a.
- No third-party MATLAB package needs to be installed: the baseline routines and
  their original licenses are included in `third_party/`.

## Quick start

Download the repository, open its root folder as MATLAB's current folder, and
close any existing parallel pool explicitly before starting the benchmark.
The runner creates one process pool and reuses it throughout a batch.

Preview the complete configured plan without loading matrices or starting workers:

```matlab
plan = run_kemeny('PlanOnly', true);
```

Run a small example with two workers and all five methods:

```matlab
results = run_kemeny('Tests', {
    'minnesota', 1, 1, 1, [4, 1], [4, 1]
}, 'Workers', 2);
```

This requires at least three available CPUs (one client and two workers).
For this example `k=4` means four block tasks, which can execute in waves on
two workers. A partition count is not the same thing as the worker count.

Run the full configured experiment only when sufficient time and memory are available:

```matlab
[results, summary] = run_kemeny;
```

Each invocation writes to a fresh directory under `runs/`. To choose a directory,
pass `'OutputDir', 'runs/my-experiment'`; the runner refuses to overwrite existing
result files. Existing results are not automatically resumed after interruption.

## Configure datasets, partitions, and repetitions

Edit `TEST_CASES` at the start of `run_kemeny.m`. Each row has six columns:

```matlab
TEST_CASES = {
    'cage10',    3, 4, 4, [[4, 1], [6, 2], [8, 2]], [[8, 2], [4, 3]];
    'minnesota', 3, 4, 4, [[4, 1], [6, 2], [8, 2]], [[8, 2], [4, 3]];
};
```

| Column | Meaning |
|---|---|
| 1 | Dataset basename in `matrix/`; `.mat` is optional. Directory paths are not accepted. |
| 2 | Number of measured Direct calls. |
| 3 | Number of measured Bini-Rec calls. |
| 4 | Number of measured Bini-DirRec calls. |
| 5 | Du groups, each specified as `[k, measured repeats]`. |
| 6 | New-full-rhs groups, each specified as `[k, measured repeats]`. |

MATLAB concatenates `[[4,1],[6,2]]` into `[4 1 6 2]`; the runner explicitly
parses adjacent numbers as pairs. The equivalent two-column form
`[4,1; 6,2]` and the single-group form `[4,1]` are also accepted.
Groups run in their listed order. Repeated values of `k` remain separate groups
identified by `PlanIndex`. A repeat count of zero skips that group and its warm-up.

For normal HPC use, leave `WORKERS = []`. The runner chooses
`min(allocated CPUs - 1, largest enabled k)` workers.
`SC_BATCH_SIZE = Inf` and `FEEDBACK_BATCH_SIZE = 64` are the reference settings.
The former controls stochastic-complement RHS batches; the latter controls
additional feedback RHS batches. Both can change runtime and temporary storage.

## Experimental protocol

- Every enabled method/partition group receives one full, untimed-for-summary
  warm-up call before its measured repetitions.
- All five methods use caller-side `tic`/`toc` around one complete algorithm call.
  Matrix loading, graph preprocessing, `dissect` ordering, stationary distribution,
  pool setup, warm-up, and result writing are excluded from the measured `Seconds`.
- One process pool is reused across all groups and datasets. Each worker has one
  computational thread. The Du/New client uses one thread; Direct and the Bini
  methods may use up to the allocated CPU count. At most `min(k, workers)` workers
  execute the block tasks concurrently.
- `DirectMean` is the mean of successful finite **measured** Direct values in the
  same configuration row. Warm-up values are excluded.
- `RelativeError = abs(Kemeny - DirectMean) / max(abs(DirectMean), eps)`.
  Without a successful measured Direct value, the reference and error are `NaN`.
  Direct is a double-precision numerical reference, not an exact answer.
- Warm-up failures stop the run. Measured failures are recorded, and the overall
  run reports failure even if later methods complete.

For `Problem.A` inputs, the runner uses the unweighted adjacency pattern,
symmetrizes directed input, and retains one largest connected component.
For files containing a raw `A`, it preserves weights, diagonal, and direction and
requires irreducibility. It then applies `dissect`, forms the row-stochastic
transition matrix, and computes its stationary distribution. Dataset files
themselves are not modified. Thus the reported dimension may differ from the
original matrix dimension. See [dataset provenance](matrix/README.md).

The New kernel uses full RHS triangular solves, updates each stochastic complement
by batches, and reuses factorizations and local inverse matrices for the feedback
correction. Performance differences reflect the complete implementations.
They should not be attributed solely to the gamma formula.

## Output files

| File | Content |
|---|---|
| `plan.csv` | Expanded algorithm/partition groups and requested counts. |
| `results-partial.tsv` | Measured results checkpointed after each call. |
| `warmup-results.tsv` | Separate warm-up results. |
| `results.csv`, `results.mat` | Final measured results. |
| `summary.csv` | Per-group counts, mean/median/min/max/std time, and maximum relative error. |
| `timing-matlab.tsv` | Stage timings and client-only CPU time. |
| `configuration.mat`, `environment.mat` | Configuration and runtime metadata. |

`TOTAL` stages include their nested child stages, so do not sum both.
`ClientCPUSeconds` excludes worker processes and cannot establish total parallel
CPU utilization. A single measured repetition is an observation, not a stability
estimate. Final summaries may be absent after a scheduler kill; completed calls
are retained in the checkpoint files.

## CSU HPC submission

`submit_free.sh` is an example for CSU's `free` account, `freecpuQ` partition, and
`normal` QOS. It requests one node, one task, 15 CPUs, 100 GiB of memory, and 24 hours.
It uses the directory from which `sbatch` is called; there is no personal home path.
MATLAB is loaded with `module load MATLAB/R2021b`.

From the repository root on the cluster:

```bash
sbatch --test-only submit_free.sh
sbatch submit_free.sh
```

Check that the account, partition, memory, and requested time are accepted by your
scheduler before submitting. Other clusters require adjustments to the directives
and MATLAB module command. The example deliberately rejects a runtime account or
partition different from its configured CSU free settings.

The public script requests 24 hours because the recorded complete experiment took
19 hours 23 minutes. The downloaded historical submission script contained an
8-hour directive, which does not explain the actual completed duration; its
effective scheduler time limit was not recorded in the supplied evidence.
Do not treat either directive as a verified cluster maximum.

On HPC, the script snapshots the MATLAB sources and matrices into `runs/free-JOBID`
when the job starts. Keep the source configuration unchanged while a job is queued.
Fresh generated logs, MATLAB preferences, certificates, and worker caches must
remain local and are excluded from this publication tree.

## Reference experiment

[Reference results](reference_results/20372192/README.md) contain the numerical
outputs of job `20372192`, including the full-precision measured checkpoints.
The reference hardware was a CSU compute node with two Intel Xeon Gold 6248R
processors (24 cores per processor, nominal frequency 3.00 GHz).
The job allocation was 15 CPU cores and 100 GiB of memory, with 14 single-threaded
MATLAB workers. Allocated memory is not the node's total physical memory.

Before publication cleanup, the runner, both `src/` files, all vendored baseline
files, and all 25 matrix files matched the saved job snapshot byte for byte.
Publication cleanup changes only whole-line comments/blank lines in first-party
MATLAB files; executable lines and the default experiment plan are unchanged.
The Slurm wrapper is generalized as described above. The supplied tests were not
part of that HPC job; they are independent release checks.
See [provenance.json](reference_results/20372192/provenance.json) for source and data hashes.

## Small correctness checks

From the repository root:

```matlab
addpath('tests');
test_kemeny_feedback;
test_sc_batches;
```

These use independent small-matrix references and do not create a parallel pool.
To test the parallel paths as well, start an appropriate process pool explicitly
and call `test_kemeny_feedback(true)` and `test_sc_batches(true)`.
The existing pool must be closed explicitly before using the benchmark runner.

## Layout

```text
run_kemeny.m                 configuration and benchmark runner
submit_free.sh               portable-path CSU Slurm example
src/kemeny_feedback.m        New-full-rhs numerical kernel
src/kemeny_prepare_pool.m    optional standalone pool helper; not used by the runner
tests/                      independent small-matrix checks
third_party/                unmodified baselines and original licenses
matrix/                     input matrices and dataset attribution
reference_results/20372192/  recorded experiment outputs and provenance
LICENSE                     MIT license for original project code/documentation
THIRD_PARTY_NOTICES.md       baseline/data copyright and license notices
```

## License and attribution

Original project code and documentation are provided under the [MIT License](LICENSE),
with copyright attributed to the repository identity `Dragonfly-Home`.
Third-party code and datasets retain their own licenses and copyright notices;
the root MIT license does not replace them.
See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and [matrix/README.md](matrix/README.md).

When citing this implementation, use the repository URL
`https://github.com/Dragonfly-Home/Kemeny-and-Gamma` together with the release tag
or commit identifier used for the experiment. Cite the relevant baseline papers
and dataset sources separately. Add the associated manuscript's bibliographic
information when it is finalized; no unpublished title, author list, or DOI is
invented here.
