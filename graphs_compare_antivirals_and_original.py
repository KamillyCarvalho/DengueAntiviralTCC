from graph_from_csv_data import *
from antiviral_adm_time import day

base_path = os.path.dirname(os.path.abspath(__file__))

population_path = "z" # Options: "s", "i", "v", "z"

if(population_path == "s"):
    label_1 = "Monocitos susceptiveis/uL"            
    eixo_y_label = "Monócitos suscetíveis/\u03bcL"
elif(population_path == "i"):
    label_1 = "Monocitos infectados/uL"            
    eixo_y_label = "Monócitos infectados/\u03bcL"  
elif(population_path == "v"):
    label_1 = "Particulas virais/uL"            
    eixo_y_label = "Partículas virais/\u03bcL"  
elif(population_path == "z"):
    label_1 = "Linfocitos T/uL"            
    eixo_y_label = "Linfócitos T/\u03bcL"   

# Define the paths to the CSV files for the three datasets
csv1_ref = os.path.join(base_path, "results", "antiviral_original_time_"+ population_path+".csv")  # Path for the first dataset
csv2_ref = os.path.join(base_path, "results", "antiviral_rho_time_"+ population_path+".csv")   # Path for the second dataset
csv3_ref = os.path.join(base_path, "results", "antiviral_xi_time_"+ population_path+".csv")  # Path for the third dataset
csv4_ref = os.path.join(base_path, "results", "antiviral_psi_time_"+ population_path+".csv")  # Path for the third dataset

# Generate the graph using the three datasets and save it as a PNG file
graph_generator_from_4_csv(
    day,                # Day of antiviral administration
    csv1_ref,           # First dataset
    csv2_ref,           # Second dataset
    csv3_ref,           # Third dataset
    csv4_ref,           # Fourth dataset
    eixo_y_label,       # Label for the Y-axis
    "grafico-"+population_path+"-4data-day-"+str(day),  # Name of the output graph file
    "Tempo",        # Label for the X-axis
    label_1             # Label for the Y-axis data
)