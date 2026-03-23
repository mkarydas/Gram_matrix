# Gram Matrix — Effective Dimension

Given an analytical kernel function K(x, y) that defines inner products, this project:

1. **Discretizes** the domain into a finite set of points
2. **Builds** the Gram matrix G where `G[i, j] = K(x_i, x_j)`
3. **Computes** the effective dimension — the number of non-negligible eigenvalues of G

## Mathematical Background

A **Gram matrix** is a symmetric positive semi-definite matrix whose entries are pairwise inner products (or kernel evaluations) between a set of points. Formally, for points x₁, …, xₙ and a kernel K:

```
G[i, j] = K(x_i, x_j)
```

The **rank** of G (number of non-zero eigenvalues) reflects the intrinsic dimensionality of the function space spanned by the kernel. In practice, eigenvalues below a small threshold are treated as zero — this count is the **effective dimension**.

## Files

| File | Description |
|------|-------------|
| `gram_matrix.py` | Core library: `discretize`, `build_gram_matrix`, `effective_dimension` |
| `kernels.py` | Built-in kernels: linear, polynomial, Gaussian, Laplacian |
| `example.py` | Runnable demo comparing kernels on a 1D domain |

## Usage

```python
from kernels import gaussian_kernel
from gram_matrix import discretize, build_gram_matrix, effective_dimension

# 1. Choose discretization points
points = discretize(domain=(0.0, 1.0), n_points=100, method="grid")

# 2. Build the Gram matrix
K = gaussian_kernel(sigma=0.5)
G = build_gram_matrix(K, points)

# 3. Compute effective dimension
dim, eigenvalues = effective_dimension(G, threshold=1e-6, relative=True)
print(f"Effective dimension: {dim}")
```

### Custom kernel

Any callable `f(x, y) -> float` works:

```python
import numpy as np

def my_kernel(x, y):
    return np.cos(x - y)  # example: cosine kernel

G = build_gram_matrix(my_kernel, points)
dim, _ = effective_dimension(G)
```

### Multi-dimensional domain

```python
# 2D grid over [0,1]²
points = discretize(domain=((0, 1), (0, 1)), n_points=20, method="grid")
# → 400 points, shape (400, 2)
```

## Dependencies

- Python 3.x
- `numpy`

Install with:
```bash
pip install numpy
```

## Running the example

```bash
python example.py
```

Example output:
```
Kernel                  Effective dim        λ_max  λ_min (non-zero)
------------------------------------------------------------------------
Linear                             2        0.3383          2.82e-02
Polynomial d=2                     3        0.2151          9.28e-05
Polynomial d=5                     6        0.0877          1.90e-09
Gaussian σ=1.0                     1       49.8329          ...
Gaussian σ=0.1                    14        5.1632          ...
Gaussian σ=0.01                   50        0.5025          ...
```

A **smaller bandwidth σ** in the Gaussian kernel makes the kernel more local, increasing the effective dimension. A **higher polynomial degree** adds more independent directions, also increasing the effective dimension.
