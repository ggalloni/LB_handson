import camb
import numpy as np


def compute_offsets(ell, spectra):
    Nl = []
    for i in range(len(ell)):
        spec = spectra[:, i]
        negative_spectra = np.abs(spec[spec < 0])
        try:
            quantile = np.quantile(negative_spectra, 0.99)
        except IndexError:
            quantile = 0.0
        Nl.append(quantile)
    return np.array(Nl)


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
