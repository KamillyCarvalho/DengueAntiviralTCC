from project_libraries import *
from generate_graphs import * 

base_path = os.path.dirname(os.path.abspath(__file__))

population_id = "z" # Options: "s", "i", "v", "z"

if(population_id == "s"):
    population_path = "susceptible_monocytes_antiviral"
    label_1 = "Monocitos susceptiveis/uL"            
    eixo_y_label = "Monócitos suscetíveis/\u03bcL"
elif(population_id == "i"):
    population_path = "infected_monocytes_antiviral"
    label_1 = "Monocitos infectados/uL"            
    eixo_y_label = "Monócitos infectados/\u03bcL"  
elif(population_id == "v"):
    population_path = "viral_particles_antiviral"
    label_1 = "Particulas virais/uL"            
    eixo_y_label = "Partículas virais/\u03bcL"  
elif(population_id == "z"):
    population_path = "t_lymphocytes_antiviral"
    label_1 = "Linfocitos T/uL"            
    eixo_y_label = "Linfócitos T/\u03bcL"   

# Define the paths to the CSV files for the three datasets
csv1_ref = os.path.join(base_path, "results_csv", population_path+"_rho.csv")  # Path for the first dataset
csv2_ref = os.path.join(base_path, "results_csv", population_path+"_xi.csv")   # Path for the second dataset
csv3_ref = os.path.join(base_path, "results_csv", population_path+"_psi.csv")  # Path for the third dataset

# Generate the graph using the three datasets and save it as a PNG file
graph_generator_from_3_csv(
    csv1_ref,           # First dataset
    csv2_ref,           # Second dataset
    csv3_ref,           # Third dataset
    eixo_y_label,       # Label for the Y-axis
    "grafico-"+population_id+"-3data", # Name of the output graph file
    "Antiviral",        # Label for the X-axis
    label_1             # Label for the Y-axis data
)