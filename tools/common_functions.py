import camb
import numpy as np


def compute_offsets(ell, varcl, clref, fsky=1.0, iter=1):
    Nl = np.sqrt(np.abs(varcl - (2.0 / (2.0 * ell + 1) * clref**2) / fsky))
    for _ in range(iter):
        Nl = np.sqrt(
            np.abs(
                varcl - 2.0 / (2.0 * ell + 1) / fsky * (clref**2 + 2.0 * Nl * clref)
            )
        )
    return Nl * np.sqrt((2.0 * ell + 1) / 2.0)


def compute_theoretical_spectrum(lmax, r):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.32, ombh2=0.02237, omch2=0.1201, mnu=0.06, omk=0, tau=0.06)
    pars.InitPower.set_params(As=2.12e-9, ns=0.9651, r=r)
    pars.set_for_lmax(lmax=2500)

    pars.WantTensors = True
    pars.DoLensing = True

    results = camb.get_results(pars)
    res = results.get_cmb_power_spectra(
        CMB_unit="muK",
        lmax=lmax,
        raw_cl=True,
    )
    return res["total"][:, [0, 1, 2, 3]]
