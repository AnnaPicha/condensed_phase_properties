# plotting_script.py
import matplotlib.pyplot as plt

def plot_data(analyzed_data):
    """
    Plots the analyzed data.
    """
    print(f"Plotting data: {analyzed_data}")
    plt.plot([0, 1, 2], [0, analyzed_data, analyzed_data / 2])
    plt.title("Analyzed Data Plot")
    plt.show()
