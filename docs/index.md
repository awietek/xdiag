---
title: Home
---
![C++ Badge](https://img.shields.io/badge/C%2B%2B-00599C?logo=cplusplus&logoColor=fff&style=for-the-badge)
![Julia](https://img.shields.io/badge/-Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white)

[![Linux CI](https://github.com/awietek/xdiag/actions/workflows/linux.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/linux.yml)
[![Mac OSX CI](https://github.com/awietek/xdiag/actions/workflows/osx.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/osx.yml)
[![Intel MPI CI](https://github.com/awietek/xdiag/actions/workflows/intelmpi.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/intelmpi.yml)
[![Julia CI](https://github.com/awietek/XDiag.jl/actions/workflows/CI.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/XDiag.jl/actions/workflows/CI.yml)

![license](https://img.shields.io/badge/license-Apache%202.0-blue)

[Quick Start](quick_start.md){ .md-button .md-button--primary }
[Code on GitHub](https://github.com/awietek/xdiag){ .md-button .md-button--secondary }

## Overview
XDiag is a library for performing Exact Diagonalizations of
quantum many-body systems. Key features include optimized combinatorical
algorithms for navigating Hilbert spaces, iterative linear algebra algorithms,
 shared and distributed memory parallelization. It consist of two packages:
 
* The core C++ library [xdiag](https://github.com/awietek/xdiag)
* The convenient Julia wrapper library [XDiag.jl](https://github.com/awietek/XDiag.jl)


## Citation
Please support our work by citing XDiag and the implemented algorithms if it is used in your published research.

```bibtex
@article{Wietek2018,
  title = {Sublattice coding algorithm and distributed memory parallelization for large-scale exact diagonalizations of quantum many-body systems},
  author = {Wietek, Alexander and L\"auchli, Andreas M.},
  journal = {Phys. Rev. E},
  volume = {98},
  issue = {3},
  pages = {033309},
  numpages = {10},
  year = {2018},
  month = {Sep},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevE.98.033309},
  url = {https://link.aps.org/doi/10.1103/PhysRevE.98.033309}
}

```


## Gallery
<div class="grid cards" markdown>
- ![Image title](img/triangularwse2.png){ align=left }
- ![Image title](img/hb_chain_dynamical_sf.png){ align=left }
- ![Image title](img/hubbard_doublon.png){ align=left }
- ![Image title](img/j1j2_spectra.png){ align=left }
</div>
