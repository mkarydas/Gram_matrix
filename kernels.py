"""
kernels.py — Common kernel (inner product) functions for building Gram matrices.

Each kernel function takes two points x, y (scalars or NumPy arrays)
and returns a scalar K(x, y).
"""

import numpy as np


def linear_kernel(x, y):
    """K(x, y) = x · y"""
    return np.dot(x, y)


def polynomial_kernel(c=0, d=2):
    """K(x, y) = (x · y + c)^d"""
    def kernel(x, y):
        return (np.dot(x, y) + c) ** d
    kernel.__name__ = f"polynomial_kernel(c={c}, d={d})"
    return kernel


def gaussian_kernel(sigma=1.0):
    """K(x, y) = exp(-||x - y||^2 / (2 * sigma^2))"""
    def kernel(x, y):
        diff = np.asarray(x) - np.asarray(y)
        return np.exp(-np.dot(diff, diff) / (2 * sigma ** 2))
    kernel.__name__ = f"gaussian_kernel(sigma={sigma})"
    return kernel


def laplacian_kernel(sigma=1.0):
    """K(x, y) = exp(-||x - y||_1 / sigma)"""
    def kernel(x, y):
        diff = np.asarray(x) - np.asarray(y)
        return np.exp(-np.sum(np.abs(diff)) / sigma)
    kernel.__name__ = f"laplacian_kernel(sigma={sigma})"
    return kernel
