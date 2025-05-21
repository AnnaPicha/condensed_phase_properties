import condensed_phase_properties.analysis.cp_props as analysis
import pandas as pd
from importlib import resources  # Python ≥ 3.9
from openff.units import unit
import sys

def load_initial_properties(species: str):
    file_name_props = f"cp_props_{species}.csv"
    data_path_props = resources.files(f"condensed_phase_properties.data.all_property_data.{species}.thermodynamic_properties.initial_npt_run") / file_name_props
    with data_path_props.open("r", encoding="utf-8") as f:
        props = pd.read_csv(f, sep="&")
    file_name_timeseries = f"cp_props_time_series_{species}.csv"
    data_path_timeseries = resources.files(f"condensed_phase_properties.data.all_property_data.{species}.thermodynamic_properties.initial_npt_run") / file_name_timeseries
    with data_path_timeseries.open("r", encoding="utf-8") as f:
        timeseries = pd.read_csv(f, sep="&")
    return props, timeseries

def load_repetition_properties(species: str, method: str):
    file_name = f"{method}_props_per_run.csv"
    data_path = resources.files(f"condensed_phase_properties.data.all_property_data.{species}.thermodynamic_properties.5_repetition_runs") / file_name
    with data_path.open("r", encoding="utf-8") as f:
        return pd.read_csv(f, sep="&")
    
def load_msd_diffusion(species: str, method: str):
    file_name_diffusion = f"diffusion_{species}_{method}_NVT_NH.dat"
    data_path = resources.files(f"condensed_phase_properties.data.all_property_data.{species}.diffusion") / file_name_diffusion
    with data_path.open("r", encoding="utf-8") as f:
        diffusion = pd.read_csv(f, sep="&")
    file_name_msd = f"msd_{species}_{method}_NVT_NH.dat"
    data_path = resources.files(f"condensed_phase_properties.data.all_property_data.{species}.diffusion") / file_name_msd
    with data_path.open("r", encoding="utf-8") as f:
        msd = pd.read_csv(f, sep="&")
    return msd, diffusion


species = sys.argv[1]

props, timeseries = load_initial_properties(species)

print("\n")
print(f"Properties from initial NPT 1ns run for {species}")
print(props)
print("-"*50)
print("\n")

for theory in ['mm', 'ani2x', 'mace_s']:
    water_5_reps = load_repetition_properties(species, theory)
    print(f"Properties from 5 repitition NPT runs for {species} with {theory}")
    print(water_5_reps)
    print("-"*50)
    print("\n")

for theory in ['mm', 'ani2x', 'mace_s']:
    msd, diffusion = load_msd_diffusion(species, theory)
    print(f"Diffusion coefficient for {theory} {species}")
    print(diffusion.columns[0])
    print("-"*50)

print("\n")