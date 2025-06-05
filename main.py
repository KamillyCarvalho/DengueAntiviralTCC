from project_libraries import *  # Importing necessary libraries for numerical calculations and data handling
from generate_graphs import * 
from config_variables import *
from dengue_math_models import *
from graph_from_csv_data import *

antiviral_RHO = "rho" 
antiviral_XI = "xi"
antiviral_PSI = "psi"
no_antiviral = "original"

antiviral_parameters = [
    antiviral_RHO, 
    antiviral_XI,
    antiviral_PSI,
]

folder_name_for_figures = "results_figures" # Folder to save the figures
folder_name_for_csv = "results_csv"         # Folder to save the CSV files
folder_name_for_csv_sensivity = "results_csv_sensivity"         # Folder to save the CSV files
folder_name_for_csv_no_antiviral = "results_csv_no_antiviral" # Folder to save the CSV files for the original model

population_id = ["s", # susceptible monocytes
                 "i", # infected monocytes
                 "v", # viral particles
                 "z", # T lymphocytes
]

time_inicial = time.time() # Start the timer

# Starting the simulation for the antiviral behavior graph
print("\nStarting the simulation for antivirals behavior graph...\n")

antiviral_behavior_graph(response, percent_doses, data, label,folder_name_for_figures)
print("Simulation for antivirals behavior graph completed!")

# Starting the simulation for the original model
print("\nStarting the simulation for antivirals sensivity analysis...\n")

for antiviral in antiviral_parameters:
    model_1(antiviral,folder_name_for_figures,folder_name_for_csv_sensivity)
    print("Antiviral: ", antiviral, " - Completed")

print("\nGenerating graphs for antiviral sensivity analysis...\n")
for population in population_id:
    generate_graphs_sensivity(population,folder_name_for_figures,folder_name_for_csv_sensivity)
    print("Figures for population: ", population, " - Completed")

# Starting the simulation for the antiviral by time (day) administration
antiviral_days = [0.5,1,1.27,2,3,4,5,6] # construir esse vetor p dias

print("\nStarting the simulation for antiviral by time with no antiviral...\n")
model_3(no_antiviral,folder_name_for_figures,folder_name_for_csv_no_antiviral)

for day in antiviral_days:
    print("\nStarting the simulation for antiviral by time administration at \n",day,"day...\n")
    folder_name_for_csv_day = folder_name_for_csv+"_day_"+ str(day)         # Folder to save the CSV files
    for antiviral in antiviral_parameters: 
        model_2(antiviral,folder_name_for_figures,folder_name_for_csv_day,day)
        print("Antiviral: ", antiviral, "for day",day, " - Completed")

    print("\nGenerating graphs to compare antivirals by day...\n")
    for population in population_id:
        generate_graphs_to_compare_antiviral_and_original(population,folder_name_for_figures,folder_name_for_csv_day,folder_name_for_csv_no_antiviral,day)
        print("Figures for population: ", population,"for day",day, " - Completed")

final_inicial = time.time() # Finishes the timer
print("\nSimulation completed in", round((final_inicial - time_inicial), 2), "seconds.\n")