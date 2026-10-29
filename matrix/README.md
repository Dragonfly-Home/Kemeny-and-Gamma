# Benchmark matrices

This directory contains the 25 MAT files preserved from the experiment
snapshot. The published batch configuration uses 24 of them; `usroads-48.mat`
is an additional dataset and was not used in that batch. The MAT files are
unchanged, including their embedded metadata. The root MIT license does not
relicense these third-party data.

## SuiteSparse Matrix Collection

The 24 matrices listed below are from the
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
| `L-9.mat` | [AG-Monien/L-9](https://sparse.tamu.edu/AG-Monien/L-9) |
| `USpowerGrid.mat` | [Pajek/USpowerGrid](https://sparse.tamu.edu/Pajek/USpowerGrid) |
| `as-22july06.mat` | [Newman/as-22july06](https://sparse.tamu.edu/Newman/as-22july06) |
| `big.mat` | [Gaertner/big](https://sparse.tamu.edu/Gaertner/big) |
| `ca-CondMat.mat` | [SNAP/ca-CondMat](https://sparse.tamu.edu/SNAP/ca-CondMat) |
| `ca-HepPh.mat` | [SNAP/ca-HepPh](https://sparse.tamu.edu/SNAP/ca-HepPh) |
| `cage10.mat` | [vanHeukelum/cage10](https://sparse.tamu.edu/vanHeukelum/cage10) |
| `cage11.mat` | [vanHeukelum/cage11](https://sparse.tamu.edu/vanHeukelum/cage11) |
| `cond-mat-2003.mat` | [Newman/cond-mat-2003](https://sparse.tamu.edu/Newman/cond-mat-2003) |
| `cs4.mat` | [DIMACS10/cs4](https://sparse.tamu.edu/DIMACS10/cs4) |
| `de2010.mat` | [DIMACS10/de2010](https://sparse.tamu.edu/DIMACS10/de2010) |
| `gre_1107.mat` | [HB/gre_1107](https://sparse.tamu.edu/HB/gre_1107) |
| `hi2010.mat` | [DIMACS10/hi2010](https://sparse.tamu.edu/DIMACS10/hi2010) |
| `minnesota.mat` | [Gleich/minnesota](https://sparse.tamu.edu/Gleich/minnesota) |
| `mycielskian15.mat` | [Mycielski/mycielskian15](https://sparse.tamu.edu/Mycielski/mycielskian15) |
| `nopoly.mat` | [Gaertner/nopoly](https://sparse.tamu.edu/Gaertner/nopoly) |
| `pesa.mat` | [Gaertner/pesa](https://sparse.tamu.edu/Gaertner/pesa) |
| `ri2010.mat` | [DIMACS10/ri2010](https://sparse.tamu.edu/DIMACS10/ri2010) |
| `tuma1.mat` | [GHS_indef/tuma1](https://sparse.tamu.edu/GHS_indef/tuma1) |
| `usps_norm_5NN.mat` | [ML_Graph/usps_norm_5NN](https://sparse.tamu.edu/ML_Graph/usps_norm_5NN) |
| `usroads-48.mat` | [Gleich/usroads-48](https://sparse.tamu.edu/Gleich/usroads-48) |
| `vsp_c-30_data_data.mat` | [DIMACS10/vsp_c-30_data_data](https://sparse.tamu.edu/DIMACS10/vsp_c-30_data_data) |
| `wing_nodal.mat` | [DIMACS10/wing_nodal](https://sparse.tamu.edu/DIMACS10/wing_nodal) |
| `worms20_10NN.mat` | [ML_Graph/worms20_10NN](https://sparse.tamu.edu/ML_Graph/worms20_10NN) |

Recommended collection citation:

Timothy A. Davis and Yifan Hu. 2011. The University of Florida Sparse Matrix
Collection. *ACM Transactions on Mathematical Software* 38(1), Article 1.
[DOI: 10.1145/2049662.2049663](https://doi.org/10.1145/2049662.2049663).

For the website interface:

Scott P. Kolodziej et al. 2019. The SuiteSparse Matrix Collection Website
Interface. *Journal of Open Source Software* 4(35), 1244.
[DOI: 10.21105/joss.01244](https://doi.org/10.21105/joss.01244).

## Pisa.mat

`Pisa.mat` contains an array named `A`, without a SuiteSparse `Problem`
metadata structure. It is reproduced unchanged from
[Kemeny-and-Conquer, commit 86a9ca8c665a255745b1d225fbd89e882a638f3e](https://github.com/Cirdans-Home/Kemeny-and-Conquer/blob/86a9ca8c665a255745b1d225fbd89e882a638f3e/matrix/Pisa.mat).
The identical file is also distributed in
[Kemeny-and-Parallel](https://github.com/Dragonfly-Home/Kemeny-and-Parallel/blob/a8b948990eb564370de86e7deb4a774cf8701a37/matrix/Pisa.mat).

Its SHA-256 is
`55afd33a60e958dc0547490e2e0cbf65a6134b7996a51e8d8fd430122629cf27`.
The originating repository's BSD 3-Clause license is retained in
[third_party/bini/LICENSE](../third_party/bini/LICENSE), with attribution to
Copyright (c) 2023, Fabio Durastante. This file is not represented as a
SuiteSparse download and is not relicensed under the root MIT license.

See [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) for the accompanying
paper and code attributions.
