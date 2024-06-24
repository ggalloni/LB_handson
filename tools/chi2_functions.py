import numpy as np


def get_fiducial_gaussian_chi2(
    data, template, noise, inv_cov, *, marginalize_cov=True
):
    Nsim = data.shape[0]
    if marginalize_cov:
        return -2 * np.log(
            (
                1
                + np.diag((data - (template + noise)[None, :]) @ inv_cov @ (data - (template + noise)[None, :]).T)
                / (Nsim - 1)
            )
            ** (-Nsim / 2)
        )
    return np.diag((data - (template + noise)[None, :]) @ inv_cov @ (data - (template + noise)[None, :]).T)


def get_ideal_ptep_chi2(data, template, noise, fsky, ell):
    M = (data) / (template + noise)[None, :]

    X = M - np.log(M) - 1
    np.nan_to_num(X, copy=False, nan=0, posinf=0, neginf=0)

    return np.sum((2 * ell[None, :] + 1) * fsky * X, axis=1)


def get_hl_chi2(
    data,
    template,
    fiducial,
    noise,
    inv_cov,
    *,
    marginalize_cov=True
):
    Nsim = data.shape[0]
    M = (data) / (template + noise)[None, :]

    ghl = np.sign(M - 1) * np.sqrt(2 * (M - np.log(M) - 1))

    X = (fiducial + noise)[None, :] * ghl
    np.nan_to_num(X, copy=False, nan=0, posinf=0, neginf=0)

    if marginalize_cov:
        return -2 * np.log((1 + np.diag(X @ inv_cov @ X.T) / (Nsim - 1)) ** (-Nsim / 2))
    return np.diag(X @ inv_cov @ X.T)


def get_offset_hl_chi2(
    data, template, fiducial, noise, offset, inv_cov, *, marginalize_cov=True
):
    Nsim = data.shape[0]
    M = (data + offset[None, :]) / (template + noise + offset)[None, :]

    ghl = np.sign(M - 1) * np.sqrt(2 * (M - np.log(np.abs(M)) - 1))

    X = (fiducial + noise + offset)[None, :] * ghl
    np.nan_to_num(X, copy=False, nan=0, posinf=0, neginf=0)

    if marginalize_cov:
        return -2 * np.log((1 + np.diag(X @ inv_cov @ X.T) / (Nsim - 1)) ** (-Nsim / 2))
    return np.diag(X @ inv_cov @ X.T)


def get_lollipop_chi2(
    data, template, fiducial, noise, offset, inv_cov, *, marginalize_cov=True
):
    Nsim = data.shape[0]
    M = (data + offset[None, :]) / (template + noise + offset)[None, :]

    glolli = (
        np.sign(M)
        * np.sign(np.abs(M) - 1)
        * np.sqrt(2 * (np.abs(M) - np.log(np.abs(M)) - 1))
    )

    X = (fiducial + noise + offset)[None, :] * glolli
    np.nan_to_num(X, copy=False, nan=0, posinf=0, neginf=0)

    if marginalize_cov:
        return -2 * np.log((1 + np.diag(X @ inv_cov @ X.T) / (Nsim - 1)) ** (-Nsim / 2))
    return np.diag(X @ inv_cov @ X.T)
