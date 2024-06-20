# Likelihood approximations for LiteBIRD

In this repo, I will show how to analyze the effects of different likelihood
approximations on a mock dataset. The main goal is to understand how these
approximations affect the constraints on the tensor-to-scalar ratio.

## Installation

To clone this repo, run

```bash
git clone https://github.com/ggalloni/LB_handson_2024.git
```

To create a conda environment (named lb_like_env) with the required packages, run

```bash
conda env create -f environment.yml
```

To activate the environment, run

```bash
conda activate lb_like_env
```

## Content

The repo already contains some necessary files:

- `masks/`: contains some masks that will be used to mimic the cut-sky case.
- `spectra/`: contains some spectra computed with `CAMB` that will be used to generate the
            mock dataset and compute the likelihoods.
- `tools/`: contains some useful functions to compute the likelihoods and generate the mock dataset.

Then, of course, there are the notebooks:

- `setting_up_the_analysis.ipynb`: this notebook will generate the simulations and
                                    compute the angular power spectra.
- `compute_the_chi2s.ipynb`: this notebook will compute the $\chi^2$ for the different
                                likelihood approximations.
- `comparing_approximations.ipynb`: this notebook will plot the results for all the
                                    approximations in order to compare them.

## Preliminary steps

Before running the notebooks, you need to install the xQML code (see [here](https://gitlab.in2p3.fr/xQML/xQML)). To do so, you need to
extract the `xQML-master.zip`, cloned from the GitLab repo. Then, after going inside the unzipped folder, you
need to run

```bash
pip install -e .
```
