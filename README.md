# Condensed phase properties

<!-- 
[![PyPI version](https://badge.fury.io/py/condensed-phase-properties.svg)](https://badge.fury.io/py/condensed-phase-properties)
![versions](https://img.shields.io/pypi/pyversions/condensed-phase-properties.svg)
-->

[![GitHub license](https://img.shields.io/github/license/AnnaPicha/condensed_phase_properties.svg)](https://github.com/AnnaPicha/condensed_phase_properties/blob/main/LICENSE)



[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


This repo condensed_phase_properties contains the functions and exemplary input data for water that were used to compute thermodynamic properties. In the data folder, all published data is available (thermodynamic data and diffusion constants).


- Free software: MIT
- Documentation: https://github.com/AnnaPicha/condensed_phase_properties/.


## Features

The functions computing condensed phase properties use openff-toolkit, which has not been released for python 3.13 yet. So create an environment with python==3.12. Then, clone and install this package with
-  pip install .
- conda install openff-toolkit -c conda-forge

All functions used to compute thermodynamic properties can be found in the analysis/cp_props.py script. In the data folder, there is some examplary data from a NPT water simulation using TIP3 water. In the data/all_property_data folder, all raw data (thermodynamic properties and diffusion constants) can be found.

## Examplary scripts

The python script `compute_properties.py` uses the water report in the data folder and the property functions in the `cp_props.py` script to compute and print thermodynamic properties. The python script `get_available_properties.py` collects and prints the raw data (stored in the data folder) for a given species. To execute those scripts (once installed this repository), run

- `python compute_properties.py`

or

- `python get_available_properties.py water`

water can be replaced by methanol, acetone, nma, benzene or hexane.

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [`mgancita/cookiecutter-pypackage`](https://mgancita.github.io/cookiecutter-pypackage/) project template.

