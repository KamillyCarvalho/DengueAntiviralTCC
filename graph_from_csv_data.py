from project_libraries import *
from generate_graphs import * 

def graph_3data_generator(csv1, csv2, csv3, eixo_y_label, nome,eixo_x,eixo_y):
    
    data1 = pd.read_csv(csv1)  # Reads data from CSV files
    data2 = pd.read_csv(csv2)
    data3 = pd.read_csv(csv3)

    x1, y1 = data1[eixo_x], data1[eixo_y]  # Assumes the CSV files have 'x' and 'y' columns
    x2, y2 = data2[eixo_x], data2[eixo_y]
    x3, y3 = data3[eixo_x], data3[eixo_y]

    fig, ax = plt.subplots()  # Creates the figure and axes
    ax.set_facecolor("#FFFFFF")
    eixo_x_label = "Antiviral (em %)"
    plt.subplots_adjust(right=0.97, top=0.97)
    plt.xlabel(eixo_x_label, fontsize=12, labelpad=4)
    plt.ylabel(eixo_y_label, fontsize=12, labelpad=6)

    plt.plot(x1, y1, '#00FF00', label='Antiviral 1', linewidth=1.9)  # Plots the three curves
    plt.plot(x2, y2, '#FFA500', label='Antiviral 2', linewidth=1.9)
    plt.plot(x3, y3, '#0000FF', label='Antiviral 3', linewidth=1.9)

    plt.grid(color='gray', linestyle='-', linewidth=0.15) # Adds grid and legend
    plt.legend()

    final_name = nome + '.png'      # Saves the graph
    save_file(final_name)           # Saves the graph using the save_file function
    return 0

base_path = os.path.dirname(os.path.abspath(__file__))

population_path = "viral_particles_antiviral" # "susceptible_monocytes_antiviral", "infected_monocytes_antiviral", "viral_particles_antiviral", "t_lymphocytes_antiviral"
label_1 = "Particulas virais/uL"            # "Monocitos susceptiveis/uL", "Monocitos infectados/uL", "Particulas virais/uL","Linfocitos T/uL"
eixo_y_label = "Partículas virais/\u03bcL"   # "Monócitos susceptíveis/\u03bcL", "Monócitos infectados/\u03bcL", "Partículas virais/\u03bcL","Linfócitos T/\u03bcL"

csv1_ref = os.path.join(base_path, "results", population_path+"_rho.csv")
csv2_ref = os.path.join(base_path, "results", population_path+"_xi.csv")
csv3_ref = os.path.join(base_path, "results", population_path+"_psi.csv")

graph_3data_generator(csv1_ref, csv2_ref, csv3_ref, eixo_y_label, 'grafico-v-3data',"Antiviral", label_1)