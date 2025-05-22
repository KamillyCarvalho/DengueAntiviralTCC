from project_libraries import *
from generate_graphs import * 

def generate_graphs_sensivity(population_id,folder_name_for_figures,folder_name_for_csv):
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

    base_path = os.path.dirname(os.path.abspath(__file__))
    # Define the paths to the CSV files for the three datasets
    csv1_ref = os.path.join(base_path, folder_name_for_csv, population_path+"_rho.csv")  # Path for the first dataset
    csv2_ref = os.path.join(base_path, folder_name_for_csv, population_path+"_xi.csv")   # Path for the second dataset
    csv3_ref = os.path.join(base_path, folder_name_for_csv, population_path+"_psi.csv")  # Path for the third dataset

    # Generate the graph using the three datasets and save it as a PNG file
    graph_generator_from_3_csv(
        csv1_ref,           # First dataset
        csv2_ref,           # Second dataset
        csv3_ref,           # Third dataset
        eixo_y_label,       # Label for the Y-axis
        "graph-"+population_id+"-sensitivity", # Name of the output graph file
        "Antiviral",        # Label for the X-axis
        label_1,             # Label for the Y-axis data
        folder_name_for_figures   # Folder to save the figure
    )

def generate_graphs_to_compare_antiviral_and_original(population_id,folder_name_for_figures,folder_name_for_csv,folder_name_for_csv_no_antiviral,day):
    if(population_id== "s"):
        label_1 = "Monocitos susceptiveis/uL"            
        eixo_y_label = "Monócitos suscetíveis/\u03bcL"
    elif(population_id == "i"):
        label_1 = "Monocitos infectados/uL"            
        eixo_y_label = "Monócitos infectados/\u03bcL"  
    elif(population_id == "v"):
        label_1 = "Particulas virais/uL"            
        eixo_y_label = "Partículas virais/\u03bcL"  
    elif(population_id == "z"):
        label_1 = "Linfocitos T/uL"            
        eixo_y_label = "Linfócitos T/\u03bcL"   

    base_path = os.path.dirname(os.path.abspath(__file__))
    # Define the paths to the CSV files for the three datasets
    csv1_ref = os.path.join(base_path, folder_name_for_csv_no_antiviral, "antiviral_original_time_"+ population_id+".csv")  # Path for the first dataset
    csv2_ref = os.path.join(base_path, folder_name_for_csv, "antiviral_rho_time_"+ population_id+"_day"+str(day)+".csv")   # Path for the second dataset
    csv3_ref = os.path.join(base_path, folder_name_for_csv, "antiviral_xi_time_"+ population_id+"_day"+str(day)+".csv")  # Path for the third dataset
    csv4_ref = os.path.join(base_path, folder_name_for_csv, "antiviral_psi_time_"+ population_id+"_day"+str(day)+".csv")  # Path for the third dataset

    # Generate the graph using the three datasets and save it as a PNG file
    graph_generator_from_4_csv(
        day,                # Day of antiviral administration
        csv1_ref,           # First dataset
        csv2_ref,           # Second dataset
        csv3_ref,           # Third dataset
        csv4_ref,           # Fourth dataset
        eixo_y_label,       # Label for the Y-axis
        "graph-"+population_id+"-compare-day-"+str(day),  # Name of the output graph file
        "Tempo",        # Label for the X-axis
        label_1,             # Label for the Y-axis data
        folder_name_for_figures # Folder to save the figure
    )