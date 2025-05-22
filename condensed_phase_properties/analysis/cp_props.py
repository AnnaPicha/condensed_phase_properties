"""Analysis package for condensed phase properties."""
import numpy as np
from numpy.typing import NDArray
from openff.units import unit
from pint import Quantity
from typing import cast
from collections import namedtuple


Constants = namedtuple('Constants', ['GAS_CONSTANT', 'BOLTZMANN_CONSTANT'])

CONSTANTS = Constants(
    GAS_CONSTANT=8.31446261815324 * unit.joule / unit.mole / unit.kelvin,
    BOLTZMANN_CONSTANT=1.380649e-23 * unit.joule / unit.kelvin,
)


def calc_heat_capacity_units(
    total_energy: NDArray[np.float64],
    number_particles: int,
    temp: float,
    molar_mass: float,
    printing: bool
) -> Quantity:
    """
    Compute the heat capacity.

    C_p = Variance(Total energy) / (number_particles * temp^2 * gas constant)
    the return value is in units cal/mole/kelvin
    for a correct unit transformation, the molar mass is required
    """
    pot_var = total_energy.var() * (unit.kilojoule / unit.mole) ** 2
    temp = temp * unit.kelvin

    val = cast(Quantity,
               pot_var
               / number_particles
               / temp**2
               / CONSTANTS.GAS_CONSTANT
               )
    val = val.to(unit.cal / unit.mole / unit.kelvin) / molar_mass
    if printing:
        print("heat capacity: ", val)
    return val


def calc_thermal_expansion(
        total_energy: NDArray[np.float64],
        volume: NDArray[np.float64],
        temp: float,
        printing: bool
) -> Quantity:
    """
    Compute the coefficient of thermal expansion.

    alpha = Cov(energy, vol) / (vol * temp^2 * gas constant)
    the return value is in units 1/Kelvin
    """
    cov_en_vol = cast(Quantity,
                      np.cov(total_energy, volume)[0][1]
                      * (unit.nanometer**3)
                      * unit.kJ / unit.mole
                      )

    T = temp * unit.kelvin
    volume = volume.mean() * (unit.nanometer**3)

    alpha = cov_en_vol / CONSTANTS.GAS_CONSTANT / T**2 / volume
    alpha_shift = cast(Quantity, alpha.to(1 / unit.kelvin))
    if printing:
        print("thermal expansion: ", alpha_shift)
    return alpha_shift


def calc_isothermal_compressibility(
        volume: NDArray[np.float64],
        temp: float,
        printing: bool
) -> Quantity:
    """
    Compute the isothermal compressibility.

    kappa = Variance(Box volume) / (k_B * temperature * volume)
    the return value is in units 1/bar
    """
    volume_var = volume.var() * (unit.nanometer**3) ** 2
    volume_mean = volume.mean() * (unit.nanometer**3)
    T = temp * unit.kelvin

    val = volume_var / CONSTANTS.BOLTZMANN_CONSTANT / T / volume_mean
    val = val.to(1 / unit.bar)
    if printing:
        print("thermal expansion: ", val)
    return val


def calc_heat_of_vaporization(
    pot_energy: NDArray[np.float64],
    pot_energy_mono: NDArray[np.float64],
    temp_traj: NDArray[np.float64],
    box_count: int,
    printing: bool
) -> Quantity:
    """
    Compute the heat of vaporization.

    Delta H_vap = mean_energy_gas - mean_energy_liquid + R*temperature
    the return value is in units kJ
    """
    pot_mean = pot_energy.mean() * unit.kilojoule / unit.mole / box_count
    pot_mono_mean = pot_energy_mono.mean() * unit.kilojoule / unit.mole
    temp_mean = temp_traj.mean() * unit.kelvin
    val = pot_mono_mean - pot_mean + CONSTANTS.GAS_CONSTANT * temp_mean
    if printing:
        print("heat of vaporozation: ", val)
    return val


def my_bootstrap_hov(
        liquid_pot: float,
        mono_pot: float,
        liquid_temp: float,
        Nboot: int,
        statfun
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    liquid_pot = np.array(liquid_pot)
    liquid_temp = np.array(liquid_temp)
    mono_pot = np.array(mono_pot)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(liquid_pot), len(liquid_pot))
        sample_pot = liquid_pot[index]
        sample_temp = liquid_temp[index]
        bastatistics = (
            statfun(sample_pot,
                    mono_pot,
                    sample_temp,
                    False).magnitude
        )
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)


def my_bootstrap_hcap(
        liquid_total: float,
        box_count: int,
        liquid_temp: float,
        Nboot: int,
        statfun
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    liquid_total = np.array(liquid_total)
    liquid_temp = np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(liquid_total), len(liquid_total))
        sample_pot = liquid_total[index]
        sample_temp = liquid_temp[index]
        bastatistics = (
            statfun(sample_pot,
                    box_count,
                    sample_temp.mean(),
                    False).magnitude
        )
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)


def my_bootstrap_texp(
        liquid_total: float,
        box_vol: float,
        liquid_temp: float,
        Nboot: int,
        statfun
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    liquid_total = np.array(liquid_total)
    box_vol = np.array(box_vol)
    liquid_temp = np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(liquid_total), len(liquid_total))
        sample_pot = liquid_total[index]
        sample_box_vol = box_vol[index]
        sample_temp = liquid_temp[index]
        bastatistics = (
            statfun(sample_pot,
                    sample_box_vol,
                    sample_temp.mean(),
                    False).magnitude
        )
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)


def my_bootstrap_icomp(
        box_vol: float,
        liquid_temp: float,
        Nboot: int,
        statfun
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    box_vol = np.array(box_vol)
    liquid_temp = np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(box_vol), len(box_vol))
        sample_box_vol = box_vol[index]
        sample_temp = liquid_temp[index]
        icomp = statfun(sample_box_vol, sample_temp.mean(), False).magnitude
        resampled_stat.append(icomp)

    return np.array(resampled_stat)
