# Third-party notices

The root [MIT license](LICENSE) applies to original Kemeny-and-Gamma code and
documentation. The third-party code and matrices below retain their own
copyright notices and license terms.

## Bini baseline implementations

The following files are reproduced without modification from
[Kemeny-and-Conquer](https://github.com/Cirdans-Home/Kemeny-and-Conquer),
commit [`86a9ca8c665a255745b1d225fbd89e882a638f3e`](https://github.com/Cirdans-Home/Kemeny-and-Conquer/tree/86a9ca8c665a255745b1d225fbd89e882a638f3e):

| Local file | Upstream file |
| --- | --- |
| `third_party/bini/kemenydirect.m` | `src/kemenydirect.m` |
| `third_party/bini/recursivekemeny.m` | `src/recursivekemeny.m` |
| `third_party/bini/recursivekemenydirect.m` | `src/recursivekemenydirect.m` |

Copyright (c) 2023, Fabio Durastante. These files are distributed under the
BSD 3-Clause License, reproduced in [third_party/bini/LICENSE](third_party/bini/LICENSE).
The original implementations, including their solver tolerances and recursive
calls, are retained for the benchmark comparison.

The upstream repository accompanies:

Dario Andrea Bini, Fabio Durastante, Sooyeong Kim, and Beatrice Meini,
"On Kemeny's constant and stochastic complement,"
*Linear Algebra and its Applications* 703 (2024), 137–162.
[DOI: 10.1016/j.laa.2024.09.001](https://doi.org/10.1016/j.laa.2024.09.001).

## Du parallel implementation

`third_party/du/kemeny_k_parallel.m` is reproduced without modification from
`src/kemeny_k_parallel.m` in
[Kemeny-and-Parallel](https://github.com/Dragonfly-Home/Kemeny-and-Parallel),
commit [`a8b948990eb564370de86e7deb4a774cf8701a37`](https://github.com/Dragonfly-Home/Kemeny-and-Parallel/tree/a8b948990eb564370de86e7deb4a774cf8701a37).

Copyright (c) 2026 Shuaibing Du. The MIT License for this implementation is
reproduced in [third_party/du/LICENSE](third_party/du/LICENSE).

The upstream README identifies the accompanying paper as
"On a Decomposition Formula for Kemeny's Constant via Stochastic Complement."

### Direct solver used by Du

`third_party/du/private/kemenydirect.m` is reproduced without modification
from `src/kemenydirect.m` in the same Kemeny-and-Parallel commit. That file
originates from Kemeny-and-Conquer and retains its upstream attribution:
Copyright (c) 2023, Fabio Durastante, BSD 3-Clause License.
The complete license is reproduced in
[third_party/du/LICENSE-Kemeny-and-Conquer](third_party/du/LICENSE-Kemeny-and-Conquer).

The private directory prevents this helper from shadowing the public Bini
baseline of the same name. Moving the file into that directory does not change
its source content or license.

## Matrix data

The matrices are not licensed under this repository's MIT license.
The [matrix provenance table](matrix/README.md) identifies each source.

All 11 included matrices are from the SuiteSparse Matrix Collection and are supplied under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/),
as stated in the collection's
[license and attribution guidance](https://sparse.tamu.edu/about#license).
Their original MAT files and embedded metadata are retained. Cite the collection
and preserve any additional matrix-specific citations in that metadata.

The runner preprocesses data in memory when constructing a transition matrix;
it does not overwrite the distributed MAT files.

## MATLAB

MATLAB and Parallel Computing Toolbox are required separately and are not
distributed with this repository. Their names are used only to identify the
software required to run the experiments.
