# busecke_abernathey_2019

[![DOI](https://zenodo.org/badge/155008528.svg)](https://zenodo.org/badge/latestdoi/155008528)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jbusecke/busecke_abernathey_2019_sciadv/master)

This repository contains the analysis code used in [Busecke and Abernathey 2019, Science Advances](http://advances.sciencemag.org/content/5/1/eaav5014).

The data used in this paper can be found [here](https://figshare.com/articles/Global_Ocean_Surface_Diffusivities_derived_from_Altimetry/4928693)

We strive to make our research open and reproducible. This repo enables the user to recalculate the SMLT estimate of surface diffusivity and change parameters. All data that is needed is provided.

Some of the data (e.g. the model runs described in this paper and processed AVISO velocities) are provided in a preprocessed form due to the large data volume.

## Run and reproduce the results

1. Clone this repository

2. Install the required python modules using [conda](https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).
The provided `environment.yml` file provides the necessary modules with appropriate version numbers.

If issues arise refer to `requirements_full.txt` for a detailed output of the conda environment used at the time of publication.
(We provide this additional information due to a bug that prevented `conda env export` to export pip installed modules and an installation of xarray from source (10.9))

3. Download all files from the [figshare repository](https://figshare.com/articles/Global_Ocean_Surface_Diffusivities_derived_from_Altimetry/4928693) into a folder named `data` in the repository root.

4. Execute and modify the notebook.
If the K_mix estimates should be regenerated (`recompute = True`), it might be necessary to delete/rename the files `K_mix_corrected.nc` and `K_mix_uncorrected.nc`
