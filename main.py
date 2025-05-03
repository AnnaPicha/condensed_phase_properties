#main.py
from condensed_phase_properties.simulation.simulation_script import run_simulation
from condensed_phase_properties.analysis.analysis_script import run_analysis
from condensed_phase_properties.plotting.plotting_script import plot_data

if __name__ == "__main__":
    # Step 1: Run simulation
    input_data_path = "condensed_phase_properties/simulation/input_data/sample_data.dat"
    simulated_data = run_simulation(input_data_path)

    # Step 2: Run analysis on simulated data
    analyzed_data = run_analysis(simulated_data)

    # Step 3: Plot the analyzed data
    plot_data(analyzed_data)
