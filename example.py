"""
example.py — Demonstrates building a Gram matrix and computing the effective dimension
             for several kernels over a 1D domain [0, 1].
"""

import numpy as np
from kernels import gaussian_kernel, polynomial_kernel, linear_kernel
from gram_matrix import discretize, build_gram_matrix, effective_dimension

N = 50  # number of discretization points
domain = (0.0, 1.0)
threshold = 1e-6   # relative threshold: ignore eigenvalues < 1e-6 * lambda_max

points = discretize(domain, n_points=N, method="grid")

kernels = [
    ("Linear",            linear_kernel),
    ("Polynomial d=2",    polynomial_kernel(c=1, d=2)),
    ("Polynomial d=5",    polynomial_kernel(c=1, d=5)),
    ("Gaussian σ=1.0",    gaussian_kernel(sigma=1.0)),
    ("Gaussian σ=0.1",    gaussian_kernel(sigma=0.1)),
    ("Gaussian σ=0.01",   gaussian_kernel(sigma=0.01)),
]

print(f"{'Kernel':<22}  {'Effective dim':>14}  {'λ_max':>12}  {'λ_min (non-zero)':>18}")
print("-" * 72)

for name, kernel in kernels:
    G = build_gram_matrix(kernel, points)
    dim, eigs = effective_dimension(G, threshold=threshold, relative=True)
    nonzero_eigs = eigs[eigs > threshold * eigs[0]]
    lambda_min = nonzero_eigs[-1] if len(nonzero_eigs) > 0 else 0.0
    print(f"{name:<22}  {dim:>14}  {eigs[0]:>12.4f}  {lambda_min:>18.2e}")
