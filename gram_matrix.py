"""
gram_matrix.py — Build a Gram matrix from an analytical kernel, discretize a domain,
and compute the effective dimension via eigenvalue analysis.
"""

import numpy as np


def discretize(domain, n_points, method="grid"):
    """
    Generate discretization points over a domain.

    Parameters
    ----------
    domain : tuple
        (a, b) for 1D, or ((a1, b1), (a2, b2), ...) for multi-dimensional.
    n_points : int
        Number of points per dimension (grid) or total points (random).
    method : str
        'grid'   — uniform grid (Cartesian product for multi-D)
        'random' — uniform random samples

    Returns
    -------
    points : np.ndarray, shape (N, d) or (N,) for 1D
    """
    # Normalize domain to list of (a, b) pairs
    if isinstance(domain[0], (int, float)):
        dims = [domain]
    else:
        dims = list(domain)

    d = len(dims)

    if method == "grid":
        axes = [np.linspace(a, b, n_points) for a, b in dims]
        grids = np.meshgrid(*axes, indexing="ij")
        points = np.stack([g.ravel() for g in grids], axis=-1)
    elif method == "random":
        points = np.column_stack([
            np.random.uniform(a, b, n_points) for a, b in dims
        ])
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'grid' or 'random'.")

    # Return 1D array for scalar domains
    if d == 1:
        return points[:, 0]
    return points


def build_gram_matrix(kernel, points):
    """
    Build the Gram (kernel) matrix K where K[i, j] = kernel(points[i], points[j]).

    Parameters
    ----------
    kernel : callable
        A function kernel(x, y) -> float.
    points : array-like, shape (N,) or (N, d)
        Discretization points.

    Returns
    -------
    G : np.ndarray, shape (N, N)
        Symmetric Gram matrix.
    """
    points = np.asarray(points)
    n = len(points)
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            val = kernel(points[i], points[j])
            G[i, j] = val
            G[j, i] = val
    return G


def effective_dimension(G, threshold=1e-6, relative=True):
    """
    Compute the effective dimension of a Gram matrix as the number of
    eigenvalues that are considered non-zero.

    Parameters
    ----------
    G : np.ndarray, shape (N, N)
        Symmetric positive semi-definite Gram matrix.
    threshold : float
        Cutoff value for deciding whether an eigenvalue is non-zero.
    relative : bool
        If True,  cutoff = threshold * lambda_max  (relative to largest eigenvalue).
        If False, cutoff = threshold                (absolute).

    Returns
    -------
    dim : int
        Effective dimension (number of significant eigenvalues).
    eigenvalues : np.ndarray
        Sorted eigenvalues (descending) for inspection.
    """
    eigenvalues = np.linalg.eigvalsh(G)        # ascending order, real for symmetric G
    eigenvalues = np.sort(eigenvalues)[::-1]   # descending

    cutoff = threshold * eigenvalues[0] if relative else threshold
    dim = int(np.sum(eigenvalues > cutoff))
    return dim, eigenvalues
