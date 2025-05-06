from graph_from_csv_data import *

base_path = os.path.dirname(os.path.abspath(__file__))

population_path = "s" # Options: "susceptible_monocytes_antiviral", "infected_monocytes_antiviral", "viral_particles_antiviral", "t_lymphocytes_antiviral"
label_1 = "Monocitos susceptiveis/uL"            # Options: "Monocitos susceptiveis/uL", "Monocitos infectados/uL", "Particulas virais/uL","Linfocitos T/uL"
eixo_y_label = "Monócitos suscetíveis/\u03bcL"    # Options:"Monócitos suscetíveis/\u03bcL", "Monócitos infectados/\u03bcL", "Partículas virais/\u03bcL","Linfócitos T/\u03bcL"

# Define the paths to the CSV files for the three datasets
csv1_ref = os.path.join(base_path, "results", "antiviral_original_time_"+ population_path+".csv")  # Path for the first dataset
csv2_ref = os.path.join(base_path, "results", "antiviral_xi_time_"+ population_path+".csv")   # Path for the second dataset
csv3_ref = os.path.join(base_path, "results", "antiviral_psi_time_"+ population_path+".csv")  # Path for the third dataset

# Generate the graph using the three datasets and save it as a PNG file
graph_3data_generator(
    csv1_ref,           # First dataset
    csv2_ref,           # Second dataset
    csv3_ref,           # Third dataset
    eixo_y_label,       # Label for the Y-axis
    'grafico-s-4data',  # Name of the output graph file
    "Tempo",        # Label for the X-axis
    label_1             # Label for the Y-axis data
)