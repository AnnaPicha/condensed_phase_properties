# analysis_script.py
from .cp_props import process_data

def run_analysis(simulated_data):
    """
    Analyzes the simulated data.
    """
    analyzed_data = process_data(simulated_data)
    print(f"Analyzed data: {analyzed_data}")
    return analyzed_data
