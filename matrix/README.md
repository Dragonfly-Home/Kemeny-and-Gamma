# Benchmark matrices

This directory contains the **11 MAT files** selected by the current default
configuration in `run_kemeny.m`. These are an unchanged subset of the matrix
files in the historical experiment snapshot, including their embedded metadata.
The root MIT license does not relicense these third-party data.

## SuiteSparse Matrix Collection

The 11 matrices listed below are from the
[SuiteSparse Matrix Collection](https://sparse.tamu.edu/). Their individual
pages identify the original contributors and additional citations.

The collection distributes the matrices under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Its [license guidance](https://sparse.tamu.edu/about#license) asks redistributors
to preserve matrix metadata, attribute the source, and identify modifications.
This repository preserves the MAT files; graph preprocessing and normalization
occur in memory only. Any later distribution of transformed matrices should
identify those changes and use names distinguishing them from the originals.

| File | Collection entry |
| --- | --- |
| `as-22july06.mat` | [Newman/as-22july06](https://sparse.tamu.edu/Newman/as-22july06) |
| `big.mat` | [Gaertner/big](https://sparse.tamu.edu/Gaertner/big) |
| `ca-CondMat.mat` | [SNAP/ca-CondMat](https://sparse.tamu.edu/SNAP/ca-CondMat) |
| `ca-HepPh.mat` | [SNAP/ca-HepPh](https://sparse.tamu.edu/SNAP/ca-HepPh) |
| `cond-mat-2003.mat` | [Newman/cond-mat-2003](https://sparse.tamu.edu/Newman/cond-mat-2003) |
| `cs4.mat` | [DIMACS10/cs4](https://sparse.tamu.edu/DIMACS10/cs4) |
| `de2010.mat` | [DIMACS10/de2010](https://sparse.tamu.edu/DIMACS10/de2010) |
| `gre_1107.mat` | [HB/gre_1107](https://sparse.tamu.edu/HB/gre_1107) |
| `hi2010.mat` | [DIMACS10/hi2010](https://sparse.tamu.edu/DIMACS10/hi2010) |
| `vsp_c-30_data_data.mat` | [DIMACS10/vsp_c-30_data_data](https://sparse.tamu.edu/DIMACS10/vsp_c-30_data_data) |
| `worms20_10NN.mat` | [ML_Graph/worms20_10NN](https://sparse.tamu.edu/ML_Graph/worms20_10NN) |

Recommended collection citation:

Timothy A. Davis and Yifan Hu. 2011. The University of Florida Sparse Matrix
Collection. *ACM Transactions on Mathematical Software* 38(1), Article 1.
[DOI: 10.1145/2049662.2049663](https://doi.org/10.1145/2049662.2049663).

For the website interface:

Scott P. Kolodziej et al. 2019. The SuiteSparse Matrix Collection Website
Interface. *Journal of Open Source Software* 4(35), 1244.
[DOI: 10.21105/joss.01244](https://doi.org/10.21105/joss.01244).

See [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) for the accompanying
paper and code attributions.
