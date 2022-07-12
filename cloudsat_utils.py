# -*- coding: utf-8 -*-
"""Auxiliary functions for sattools module."""
import numpy as np


def cc_interp2d(data, X, Z, x1, x2, nx, z1, z2, nz, use_numba=True):
    res = _interp2d(data, X, Z, x1, x2, nx, z1, z2, nz)
    return res


def _interp2d(data, X, Z, x1, x2, nx, z1, z2, nz):
    """
    Interpolate 2D data with coordinates given by 1D and 2D arrays.

    data is a two-dimensional array of data to be interpolated.
    X and Z are one- and two-dimensional arrays, giving coordinates
    of data points along the first and second axis, respectively

    data, X and Z are expected to be C-contiguous float32 numpy arrays
    with no mask and no transformation (such as transposition) applied.
    """

    xs = (x2 - x1) / nx
    zs = (z2 - z1) / nz
    w = data.shape[0]
    h = data.shape[1]

    out = np.zeros((nx, nz), dtype=np.float32)
    q = np.zeros((nx, nz), dtype=np.int32)

    for i in range(w):
        n1 = ((X[i - 1] + X[i]) / 2 - x1) / xs if i - 1 >= 0 else -1
        n2 = ((X[i + 1] + X[i]) / 2 - x1) / xs if i + 1 < w else nx
        if n2 - n1 < 1:
            n1 = n2 = (X[i] - x1) / xs

        for j in range(h):
            m1 = ((Z[i, j - 1] + Z[i, j]) / 2 - z1) / zs if j - 1 >= 0 else -1
            m2 = ((Z[i, j + 1] + Z[i, j]) / 2 - z1) / zs if j + 1 < h else nz
            if m2 - m1 < 1:
                m1 = m2 = (Z[i, j] - z1) / zs

            for n in range(int(n1 + 0.5), int(n2 + 0.5 + 1)):
                for m in range(int(m1 + 0.5), int(m2 + 0.5 + 1)):
                    if n < 0 or n >= nx or m < 0 or m >= nz:
                        continue
                    if np.isnan(data[i, j]):
                        continue
                    out[n, m] += data[i, j]
                    q[n, m] += 1

    for n in range(nx):
        for m in range(nz):
            if q[n, m] == 0:
                out[n, m] = np.nan
            else:
                out[n, m] /= q[n, m]
    return out


def _interp2d_ecmwf(data, X, Z, x1, x2, nx, z1, z2, nz):
    """
    Interpolate 2D data with coordinates given by 1D and 2D arrays.

    data is a two-dimensional array of data to be interpolated.
    X and Z are one-dimensional arrays, giving coordinates
    of data points along the first axis

    data, X and Z are expected to be C-contiguous float32 numpy arrays
    with no mask and no transformation (such as transposition) applied.
    """
    xs = (x2 - x1) / nx
    zs = (z2 - z1) / nz
    w = data.shape[0]
    h = data.shape[1]

    out = np.zeros((nx, nz), dtype=np.float32)
    q = np.zeros((nx, nz), dtype=np.int32)

    for i in range(w):
        n1 = ((X[i - 1] + X[i]) / 2 - x1) / xs if i - 1 >= 0 else -1
        n2 = ((X[i + 1] + X[i]) / 2 - x1) / xs if i + 1 < w else nx
        if n2 - n1 < 1:
            n1 = n2 = (X[i] - x1) / xs

        for j in range(h):
            m1 = ((Z[j - 1] + Z[j]) / 2 - z1) / zs if j - 1 >= 0 else -1
            m2 = ((Z[j + 1] + Z[j]) / 2 - z1) / zs if j + 1 < h else nz
            if m2 - m1 < 1:
                m1 = m2 = (Z[j] - z1) / zs

            for n in range(int(n1 + 0.5), int(n2 + 0.5 + 1)):
                for m in range(int(m1 + 0.5), int(m2 + 0.5 + 1)):
                    if n < 0 or n >= nx or m < 0 or m >= nz:
                        continue
                    if np.isnan(data[i, j]):
                        continue
                    out[n, m] += data[i, j]
                    q[n, m] += 1

    for n in range(nx):
        for m in range(nz):
            if q[n, m] == 0:
                out[n, m] = np.nan
            else:
                out[n, m] /= q[n, m]
    return out