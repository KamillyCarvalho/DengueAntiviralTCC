from project_libraries import *
from generate_graphs import * 

def graph_3data_generator(csv1, csv2, csv3, eixo_y_label, nome,eixo_x,eixo_y):
    
    data1 = pd.read_csv(csv1)  # Reads data from CSV files
    data2 = pd.read_csv(csv2)
    data3 = pd.read_csv(csv3)

    x1, y1 = data1[eixo_x], data1[eixo_y]  # Assumes the CSV files have 'x' and 'y' columns
    x2, y2 = data2[eixo_x], data2[eixo_y]
    x3, y3 = data3[eixo_x], data3[eixo_y]

    fig, ax = plt.subplots(figsize=(8, 6))  # Creates the figure and axes
    ax.set_facecolor("#FFFFFF")
    eixo_x_label = "Antiviral (em %)"
    plt.subplots_adjust(right=0.98, top=0.98)
    plt.xlabel(eixo_x_label, fontsize=20, labelpad=8)
    plt.ylabel(eixo_y_label, fontsize=20, labelpad=6)

    # Aumenta o tamanho da fonte dos valores nos eixos
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    plt.plot(x1, y1, '#00FF00', label='Antiviral 1', linewidth=3)  # Plots the three curves
    plt.plot(x2, y2, '#FFA500', label='Antiviral 2', linewidth=3)
    plt.plot(x3, y3, '#0000FF', label='Antiviral 3', linewidth=3)

    plt.grid(color='gray', linestyle='--', linewidth=0.3) # Adds grid and legend
    plt.legend(fontsize=18)

    final_name = nome + '.png'      # Saves the graph
    save_file(final_name)           # Saves the graph using the save_file function
    return 0

base_path = os.path.dirname(os.path.abspath(__file__))

population_path = "t_lymphocytes_antiviral" # Options: "susceptible_monocytes_antiviral", "infected_monocytes_antiviral", "viral_particles_antiviral", "t_lymphocytes_antiviral"
label_1 = "Linfocitos T/uL"             # Options: "Monocitos susceptiveis/uL", "Monocitos infectados/uL", "Particulas virais/uL","Linfocitos T/uL"
eixo_y_label = "Linfócitos T/\u03bcL"    # Options:"Monócitos susceptíveis/\u03bcL", "Monócitos infectados/\u03bcL", "Partículas virais/\u03bcL","Linfócitos T/\u03bcL"

# Define the paths to the CSV files for the three datasets
csv1_ref = os.path.join(base_path, "results", population_path+"_rho.csv")  # Path for the first dataset
csv2_ref = os.path.join(base_path, "results", population_path+"_xi.csv")   # Path for the second dataset
csv3_ref = os.path.join(base_path, "results", population_path+"_psi.csv")  # Path for the third dataset

# Generate the graph using the three datasets and save it as a PNG file
graph_3data_generator(
    csv1_ref,           # First dataset
    csv2_ref,           # Second dataset
    csv3_ref,           # Third dataset
    eixo_y_label,       # Label for the Y-axis
    'grafico-z-3data',  # Name of the output graph file
    "Antiviral",        # Label for the X-axis
    label_1             # Label for the Y-axis data
)