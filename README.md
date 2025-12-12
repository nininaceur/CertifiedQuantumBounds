# CertifiedQuantumBounds.jl

**CertifiedQuantumBounds.jl** is a Julia package for turning *numerical* SDP certificates arising from non-commutative polynomial optimization (NPO)—in particular NPA/SOHS-type relaxations—into **exact rational, certifiable bounds** via *rounding + Frobenius-optimal projection + PSD lifting*.

The package offers certified bounds on arbitrary non-commutative, constraint polynomial optimization problems (dense and sparse), as well as on Heisenberg chain groundstate energies and -correlations over symmetry-adapted bases.

---

## Installation
Install via:
Pkg.add(url="https://github.com/nininaceur/CertifiedQuantumBounds.git")

For arbitrary non-commutative optimization problems, NCTSSOS (https://github.com/wangjie212/NCTSSOS) is used to compute lower bounds, which are then converted into certified bounds according to the *Round + Project + Lift* procedure described in the paper. For quantum many-body problems, QMBCertify (https://github.com/wangjie212/QMBCertify) is used to compute the numerical bounds.

## Example usage


### Maximal violation of Bell inequalities

```julia
include("benchmarks/bell_test.jl")
```

### Certified bounds on ground state energies of the Heisenberg $J_1-J_2$ chain:

```julia
include("benchmarks/hs_energy.jl")
```

### Certified bounds on ground state observables of the Heisenberg $J_1-J_2$ chain:

```julia
include("benchmarks/hs_correlations.jl")
```


