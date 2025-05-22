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



skip_size = 0.090909
box_count=572
molar_mass = 18.015 * unit.gram / unit.mole

theories = ['mm', 'ani2x', 'mace_s', 'mace_m']

print("\nCOMPUTING CONDENSED PHASE PROPERTIES\n")
for theory in theories:
    liquid = get_liquid_traj(theory, 1, 11)
    gas = get_gas_traj(theory)

    skip_part_gas = int(round(gas["Potential Energy (kJ/mole)"].count()*skip_size,0))
    gas_cut = gas[skip_part_gas-1:-1] # skip the first 10%

    skip_part_liquid = int(round(liquid["Potential Energy (kJ/mole)"].count()*skip_size,0))
    liquid_cut = liquid[skip_part_liquid-1:-1] # skip the first 10%
    
    heat_capacity = analysis.calc_heat_capacity_units(
        liquid_cut["Total Energy (kJ/mole)"].to_numpy(), 
        box_count, 
        liquid_cut["Temperature (K)"].mean(), 
        molar_mass, 
        False
    )
    
    thermal_expansion = analysis.calc_thermal_expansion(
        liquid_cut["Total Energy (kJ/mole)"].to_numpy(), 
        liquid_cut["Box Volume (nm^3)"].to_numpy(), 
        liquid_cut["Temperature (K)"].mean(), 
        False
    )
    
    iso_comp = analysis.calc_isothermal_compressibility(
        liquid_cut["Box Volume (nm^3)"].to_numpy(), 
        liquid_cut["Temperature (K)"].mean(), 
        False
    )
    hov = analysis.calc_heat_of_vaporization (
        liquid_cut["Potential Energy (kJ/mole)"].to_numpy(), 
        gas_cut["Potential Energy (kJ/mole)"].to_numpy(), 
        liquid_cut["Temperature (K)"].to_numpy(), 
        box_count, 
        False
    )

    density = liquid_cut["Density (g/mL)"].mean()
    
    print(f"{'Property':35} {'Model':10} {'Value':>10} {'Unit':>25}")
    print("-" * 85)
    print(f"{'Heat capacity':35} {theory:<10} {round(heat_capacity.magnitude, 2):>10} {str(heat_capacity.units):>25}")
    print(f"{'Thermal expansion (*1e2)':35} {theory:<10} {round(thermal_expansion.magnitude*1e2, 2):>10} {str(thermal_expansion.units):>25}")
    print(f"{'Isothermal compressibility (*1e4)':35} {theory:<10} {round(iso_comp.magnitude*1e4, 2):>10} {str(iso_comp.units):>25}")
    print(f"{'Heat of vaporization':35} {theory:<10} {round(hov.magnitude, 2):>10} {str(hov.units):>25}")
    print(f"{'Density':35} {theory:<10} {round(density, 2):>10} {str((unit.gram / unit.milliliter)):>25}")
    print("-" * 85)

