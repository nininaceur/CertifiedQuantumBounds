# CertifiedQuantumBounds.jl

**CertifiedQuantumBounds.jl** is a Julia package for turning *numerical* SDP certificates arising from moment relaxations for non-commutative polynomial optimization (NPO) into **exact rational, certifiable bounds** via *rounding + Frobenius-optimal projection + PSD lifting*.

The package offers certified bounds on arbitrary non-commutative, constrained polynomial optimization problems (dense and sparse), as well as on Heisenberg chain groundstate energies and -correlations over symmetry-adapted bases.

---

## Installation

### Install explicit Dependencies:
To run, install:

**NCTSSOS** (https://github.com/wangjie212/NCTSSOS)  to compute bounds to non-commutative polynomial optimization problems via (sparse) moment relaxations, which are then converted into certified bounds according to the *Round + Project + Lift* procedure described in the paper. 

```julia
Pkg.add(url="https://github.com/wangjie212/NCTSSOS")
```
**QMBCertify** (https://github.com/wangjie212/QMBCertify) to compute the numerical bounds for energies and other ground state observables of quantum many-body systems, exploiting underlying symmetries.

```julia
Pkg.add(url="https://github.com/wangjie212/QMBCertify")
```

**Mosek** (https://github.com/MOSEK/Mosek.jl)

### Install CertifiedQuantumBounds.jl:

```julia
Pkg.add(url="https://github.com/nininaceur/CertifiedQuantumBounds")
```

---

## Example usage


### Maximal violation of Bell inequalities

```julia

# CHSH, dense

@ncpolyvar X[1:2]
@ncpolyvar Y[1:2]

CHSH = -symmetrize(X[1]*Y[1] +  X[1]*Y[2]  + X[2]*Y[1]  - X[2]*Y[2])

rational_certificate(CHSH, [], [], [X;Y], 2; partition=2, constraint="unipotent", QUIET=false, QUIETTS=true, tol=10e-30)

# CHSH, sparse

newbound, oldbound, shift = rational_certificate_sparse(CHSH, Polynomial[], Polynomial[], [X;Y],1; partition  = 2, constraint = "unipotent", QUIET = false, tol = 1e-20 )
```

### Certified bounds on ground state energies of the Heisenberg $J_1-J_2$ chain:

```julia

J2   = 0.2                      
supp = [[1;4], [1;7]]             
coe  = [3/4; 3/4*J2]               
tt   = [1;1]                    

N = 20

opt, data = GSB(supp, coe, N, 2;
    QUIET=true, rdm=0, lol=N, extra=1, pso=0, lso=0, three_type=tt, Gram=true)

result = certify_qmb(data, N, coe[1], opt; tol_gram=1e-15, tol_dft=1e-12, snn=true, J2=J2)

```

### Certified bounds on ground state observables of the Heisenberg $J_1-J_2$ chain:

```julia
J2 = 0.2

res = certify_qmb_corr(
    12, # N
    3, # Relaxation order for energy
    3; # Relaxation order for correlations
    J2 = J2, # J2
    dist = 1, # Define correlator to be bounded
    extra_E = 3, # Extra monomials for energy
    extra_corr = 3, # Extra monomials for correlations
    QUIET = true,
    tol_gram = 1e-13, # Gram rounding precision
    tol_dft  = 1e-12 # DFT matrix rounding precision
)

```


