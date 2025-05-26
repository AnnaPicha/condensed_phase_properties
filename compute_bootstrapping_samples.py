import condensed_phase_properties.analysis.cp_props as analysis
import pandas as pd
from importlib import resources  # Python ≥ 3.9
from openff.units import unit

def load_csv(filename: str, method: str) -> pd.DataFrame:
    data_path = resources.files(f"condensed_phase_properties.data.water_traj.{method}") / filename
    with data_path.open("r", encoding="utf-8") as f:
        return pd.read_csv(f, sep="\t")
    

def get_liquid_traj(method: str, start: int, end: int):
    csv_complete=pd.DataFrame()
    for i in range(start,end+1): 
        csv_part = load_csv(f"{method}_tip572_{i}_NPT.csv", method)
        csv_complete = pd.concat([csv_complete, csv_part])
    return csv_complete

def get_gas_traj(method: str):
    csv_complete=load_csv(f"gas_{method}.csv", method)
    return csv_complete

def print_results(boot_hov_mean: float,
                  boot_hov_std: float,
                  boot_hcap_mean: float,
                  boot_hcap_std: float,
                  boot_texp_mean: float,
                  boot_texp_std: float,
                  boot_icomp_mean: float,
                  boot_icomp_std: float):
    print(f"{'Heat of vaporization - mean:':45} {boot_hov_mean:.6f}")
    print(f"{'Heat of vaporization - std:':45} {boot_hov_std:.6f}\n")
    
    print(f"{'Heat capacity - mean:':45} {boot_hcap_mean:.6f}")
    print(f"{'Heat capacity - std:':45} {boot_hcap_std:.6f}\n")
    
    print(f"{'Thermal expansion coeff. - mean:':45} {boot_texp_mean:.6f}")
    print(f"{'Thermal expansion coeff. - std:':45} {boot_texp_std:.6f}\n")
    
    print(f"{'Isothermal compressibility - mean:':45} {boot_icomp_mean:.6f}")
    print(f"{'Isothermal compressibility - std:':45} {boot_icomp_std:.6f}\n")
    print('-'*60)

skip_size = 0.090909
box_count=572
molar_mass = 18.015 * unit.gram / unit.mole
Nboot=100

theories = ['mm', 'ani2x', 'mace_s', 'mace_m']

print("\nCOMPUTING BOOTSTRAPPING VALUES FOR CONDENSED PHASE PROPERTIES\n")
for theory in theories:
    liquid = get_liquid_traj(theory, 1, 11)
    gas = get_gas_traj(theory)

    skip_part_gas = int(round(gas["Potential Energy (kJ/mole)"].count()*skip_size,0))
    gas_cut = gas[skip_part_gas-1:-1] # skip the first 10%

    skip_part_liquid = int(round(liquid["Potential Energy (kJ/mole)"].count()*skip_size,0))
    liquid_cut = liquid[skip_part_liquid-1:-1] # skip the first 10%


    boot_hov = analysis.my_bootstrap_hov(liquid_pot=liquid_cut["Potential Energy (kJ/mole)"],
                                               mono_pot=gas_cut["Potential Energy (kJ/mole)"],
                                               liquid_temp=liquid_cut["Temperature (K)"],
                                               box_count=box_count,
                                               Nboot=Nboot,
                                               statfun=analysis.calc_heat_of_vaporization
                                               )

    boot_hcap = analysis.my_bootstrap_hcap(liquid_total=liquid_cut["Total Energy (kJ/mole)"],
                                                 box_count=box_count,
                                                 liquid_temp=liquid_cut["Temperature (K)"],
                                                 molar_mass=molar_mass,
                                                 Nboot=Nboot,
                                                 statfun=analysis.calc_heat_capacity_units
                                                 )

    boot_texp = analysis.my_bootstrap_texp(liquid_total=liquid_cut["Total Energy (kJ/mole)"],
                                                 box_vol=liquid_cut['Box Volume (nm^3)'],
                                                 liquid_temp=liquid_cut["Temperature (K)"],
                                                 Nboot=Nboot,
                                                 statfun=analysis.calc_thermal_expansion
                                                 )

    boot_icomp= analysis.my_bootstrap_icomp(box_vol=liquid_cut['Box Volume (nm^3)'],
                                                  liquid_temp=liquid_cut["Temperature (K)"],
                                                  Nboot=Nboot,
                                                  statfun=analysis.calc_isothermal_compressibility
                                                  )
    
    print(f'Results for {theory}:')
    print('-'*60)
    print_results(boot_hov.mean(),
                  boot_hov.std(),
                  boot_hcap.mean(),
                  boot_hcap.std(),
                  boot_texp.mean(),
                  boot_texp.std(),
                  boot_icomp.mean(),
                  boot_icomp.std()
                  )
    
    