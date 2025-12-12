# CertifiedQuantumBounds.jl

**CertifiedQuantumBounds.jl** is a Julia package for turning *numerical* SDP certificates arising from non-commutative polynomial optimization (NPO)—in particular NPA/SOHS-type relaxations—into **exact rational, certifiable bounds** via *rounding + Frobenius-optimal projection + PSD lifting*.

The package offers certified bounds on arbitrary non-commutative, constraint polynomial optimization problems (dense and sparse), as well as on Heisenberg chain groundstate energies and -correlations over symmetry-adapted bases.

---

## Installation
Install via:
Pkg.add(url="https://github.com/nininaceur/CertifiedQuantumBounds.git")
