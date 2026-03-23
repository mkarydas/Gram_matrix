# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Given an analytical kernel function K(x, y), this project discretizes the domain into a finite point set, builds the Gram matrix G[i,j] = K(x_i, x_j), and computes the **effective dimension** — the number of non-negligible eigenvalues — via `numpy.linalg.eigvalsh` (Python) or `Eigenvalues` (Mathematica).

## Running the Code

```bash
# Python
python3 example.py

# Mathematica (terminal)
math -script gram_matrix.wl
```

## File Structure

| File | Purpose |
|------|---------|
| `gram_matrix.py` | Core: `discretize()`, `build_gram_matrix()`, `effective_dimension()` |
| `kernels.py` | Kernel functions: `gaussian_kernel(sigma)`, `polynomial_kernel(c,d)`, `linear_kernel`, `laplacian_kernel(sigma)` |
| `example.py` | Demo: compares 6 kernels on a 50-point grid, prints effective dim + eigenvalue info |
| `gram_matrix.wl` | Wolfram Language mirror of the Python implementation |

## Key Design Decisions

- **Vectorized kernel computation** (both Python and Mathematica): kernels are never evaluated element-wise. Python uses `numpy` array ops; Mathematica uses `DistanceMatrix` + tensor dot products — avoids n² individual function calls.
- **`effective_dimension()` returns `(dim, eigenvalues)`** so callers can inspect the full spectrum.
- **`relative=True` threshold** (default): cutoff = `threshold × λ_max`, making it scale-invariant across kernels with very different magnitudes.
- Custom kernels: pass any `f(x, y) -> float` to `build_gram_matrix()` in Python, or use `gramCustom[kernel]` in Mathematica.

## GitHub Workflow

After every change: commit with a descriptive message and push.

```bash
git add <files> && git commit -m "..." && git push origin main
```
