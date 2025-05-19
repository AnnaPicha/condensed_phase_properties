import pandas as pd
import numpy as np
from openff.units import unit



R   = 8.31446261815324 
k_B = 1.380649 *10**(-23) 

gas_constant=R*unit.joule/unit.mole/unit.kelvin
boltzmann_constant=k_B*unit.joule/unit.kelvin


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# define NPT properties -------------------------------------------------------------------------------------------------------------------------------------------------------------

def calc_heat_capacity_units(total_energy: float, number_particles : int, temp : float, molar_mass : float, printing : bool):
    '''
        This function computes the heat capacity of a homogenious liquid using data from an NPT simulation:
        C_p = Variance(Total energy) / (number_particles * temp^2 * gas constant)
        the return value is in units cal/mole/kelvin
        for a correct unit transformation, the molar mass is required
    '''

    pot_var = total_energy.var() * (unit.kilojoule / unit.mole)**2
    temp = temp*unit.kelvin

    val = pot_var/number_particles/temp**2/gas_constant 
    val = val.to(unit.cal / unit.mole / unit.kelvin) / molar_mass
    if (printing):
        print("heat capacity: ", val)
    return val

def calc_thermal_expansion(total_energy: float, volume: float, temp: float, printing : bool):
    '''
        This function computes the coefficient of thermal expansion of a homogenious liquid using data from an NPT simulation:
        alpha = Covariance(Total energy, box volume) / (box volume * temp^2 * gas constant)
        the return value is in units 1/Kelvin
    '''

    cov_en_vol = np.cov(total_energy,volume)[0][1]*(unit.nanometer**3)* unit.kJ / unit.mole
    T=temp*unit.kelvin
    volume=volume.mean()*(unit.nanometer**3)
    alpha = cov_en_vol/gas_constant/T**2/volume
    alpha_shift=alpha.to(1/unit.kelvin)
    if (printing):
        print("thermal expansion: ", alpha_shift)
    return alpha_shift

def calc_isothermal_compressibility(volume : float, temp : float, printing : bool):
    '''
        This function computes the isothermal compressibility of a homogenious liquid using data from an NPT simulation:
        kappa = Variance(Box volume) / (k_B * temperature * volume)
        the return value is in units 1/bar
    '''

    volume_var = volume.var()*(unit.nanometer**3)**2
    volume_mean = volume.mean()*(unit.nanometer**3)
    T=temp*unit.kelvin

    val = volume_var/boltzmann_constant/T/volume_mean
    val = val.to(1/unit.bar)
    if (printing):
        print("thermal expansion: ", val)
    return val

def calc_heat_of_vaporization (pot_energy: float, pot_energy_mono: float, temp_traj: float, box_count : int, printing: bool):
    '''
        This function computes the heat of vaporization of a homogenious liquid using data from an NPT simulation:
        Delta H_vap = mean_energy_gas - mean_energy_liquid + R*temperature
        the return value is in units kJ
    '''

    pot_mean = pot_energy.mean() * unit.kilojoule/unit.mole / box_count
    pot_mono_mean = pot_energy_mono.mean() * unit.kilojoule/unit.mole
    temp_mean = temp_traj.mean() * unit.kelvin
    val =  pot_mono_mean - pot_mean + gas_constant * temp_mean
    if (printing):
        print("heat of vaporozation: ", val)
    return val

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# define bootstrapping functions ----------------------------------------------------------------------------------------------------------------------------------------------------

def my_bootstrap_hov(liquid_pot : float, mono_pot : float, liquid_temp : float, Nboot : int, statfun):
    '''Calculate bootstrap statistics for a sample x'''
    liquid_pot=np.array(liquid_pot)
    liquid_temp=np.array(liquid_temp)
    mono_pot = np.array(mono_pot)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0,len(liquid_pot),len(liquid_pot))
        sample_pot = liquid_pot[index]
        sample_temp = liquid_temp[index]
        bastatistics = statfun(sample_pot, mono_pot, sample_temp, False).magnitude
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)

def my_bootstrap_hcap(liquid_total : float, box_count : int, liquid_temp : float, Nboot : int, statfun):
    '''Calculate bootstrap statistics for a sample x'''
    liquid_total=np.array(liquid_total)
    liquid_temp=np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0,len(liquid_total),len(liquid_total))
        sample_pot = liquid_total[index]
        sample_temp = liquid_temp[index]
        bastatistics = statfun(sample_pot, box_count, sample_temp.mean(), False).magnitude
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)

def my_bootstrap_texp(liquid_total : float, box_vol : float, liquid_temp : float, Nboot : int, statfun):
    '''Calculate bootstrap statistics for a sample x'''
    liquid_total=np.array(liquid_total)
    box_vol=np.array(box_vol)
    liquid_temp=np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0,len(liquid_total),len(liquid_total))
        sample_pot = liquid_total[index]
        sample_box_vol = box_vol[index]
        sample_temp = liquid_temp[index]
        bastatistics = statfun(sample_pot, sample_box_vol, sample_temp.mean(), False).magnitude
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)

def my_bootstrap_icomp(box_vol : float, liquid_temp : float, Nboot : int, statfun):
    '''Calculate bootstrap statistics for a sample x'''
    box_vol=np.array(box_vol)
    liquid_temp=np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0,len(box_vol),len(box_vol))
        sample_box_vol = box_vol[index]
        sample_temp = liquid_temp[index]
        icomp = statfun(sample_box_vol, sample_temp.mean(), False).magnitude
        resampled_stat.append(icomp)

    return np.array(resampled_stat)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


