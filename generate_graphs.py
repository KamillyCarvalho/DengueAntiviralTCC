from project_libraries import *  # Importing necessary libraries for numerical calculations and data handling

def save_file(nome,folder_name):                       # Saves the file in a specific results folder
    output_dir = folder_name        # Name of the output folder
    if not os.path.exists(output_dir):     # Checks if the folder exists
        os.makedirs(output_dir)            # Creates the folder if it does not exist
    name = os.path.join(output_dir, nome)  # Full path of the file
    plt.savefig(name, format='png')        # Saves the graph in PNG format
    return 0

def save_data_to_csv(x_data, y_data, file_name, x_label, y_label,folder_name): #Saves the data being plotted to a CSV file.
    output_dir = folder_name              # Directory to save the CSV file
    if not os.path.exists(output_dir):  # Check if the directory exists
        os.makedirs(output_dir)         # Create the directory if it doesn't exist

    file_path = os.path.join(output_dir, file_name + ".csv")  # Full path for the CSV file
   
    with open(file_path, "w") as file:        # Save the data to a CSV file
        file.write(f"{x_label},{y_label}\n")  # Write the header
        for x, y in zip(x_data, y_data):
            file.write(f"{x},{y}\n")          # Write the data rows

def graph_generator_3(eixo_x,eixo_y,t,dados,dados2,nome,folder_name,day_antiviral):
    fig, ax = plt.subplots()                         # Creates a figure and an axis
    ax.set_facecolor("#FFFFFF")                      # Sets the axis background to white
    plt.subplots_adjust(right=0.97, top=0.97)        # Adjusts the position of the subplots
    plt.xlabel(eixo_x, fontsize = 12, labelpad = 4)  # Sets the x-axis label
    plt.ylabel(eixo_y,  fontsize = 12, labelpad = 6) # Sets the y-axis label
    plt.plot(t,dados[0,:,0],'#00c2b0',label = 'Com antiviral',linewidth = '1.9')  # Plots the data with antiviral
    plt.plot(t,dados2[0,:,0],'#ff8201',label = 'Sem antiviral',linewidth = '1.9') # Plots the data without antiviral
    plt.axvline(day_antiviral, color = 'b', linestyle = '--', linewidth = 1.3,label = 'Adm. antiviral') # Adds a vertical line indicating the day of antiviral administration
    plt.tight_layout()
    leg = ax.legend()                                # Creates the legend
    plt.grid(color = 'gray', linestyle = '-', linewidth = 0.15) # Adds a grid to the graph
    final_name = nome + '[' + str(day_antiviral) + '].png'     # Final file name, including the day of antiviral administration
    save_file(final_name,folder_name)                            # Saves the graph using the save_file function
    plt.close(fig)  # Closes the figure to free up memory
    return 0

def graph_2D_generator(t,data1, data2, eixo_y_label, nome,folder_name):
    fig, ax = plt.subplots(figsize=(8, 6))                               # Creates a figure and an axis
    ax.set_facecolor("#FFFFFF")                          # Sets the axis background to white
    eixo_x_label = "Tempo (em dias)"                       # Sets the x-axis label
    
    plt.subplots_adjust(right=0.99, top=0.98, left=0.15, bottom=0.12)
    plt.xlabel(eixo_x_label, fontsize=20, labelpad=8)
    plt.ylabel(eixo_y_label, fontsize=20, labelpad=1)
    plt.plot(data1, data2, "#CC1100", label='Dados', linewidth=3)  # Plots the data
    plt.grid(color='gray', linestyle='-', linewidth=0.15)  # Adds a grid to the graph

    ax.tick_params(axis='x', labelsize=16) # Increase the font size of the axis tick labels
    ax.tick_params(axis='y', labelsize=16)
    plt.tight_layout()

    final_name = nome + '.png'                             # Final file name
    save_file(final_name,folder_name)                                  # Saves the graph using the save_file function
    plt.close(fig)  # Closes the figure to free up memory
    return 0

def graph_generator_from_3_csv(csv1, csv2, csv3, eixo_y_label, nome,eixo_x,eixo_y,folder_name):
    
    data1 = pd.read_csv(csv1)  # Reads data from CSV files
    data2 = pd.read_csv(csv2)
    data3 = pd.read_csv(csv3)

    x1, y1 = data1[eixo_x], data1[eixo_y]  # Assumes the CSV files have 'x' and 'y' columns
    x2, y2 = data2[eixo_x], data2[eixo_y]
    x3, y3 = data3[eixo_x], data3[eixo_y]

    fig, ax = plt.subplots(figsize=(8, 6))  # Creates the figure and axes
    ax.set_facecolor("#FFFFFF")
    eixo_x_label = "Eficácia antiviral (em %)"
    # plt.subplots_adjust(right=0.99, top=0.98)
    plt.subplots_adjust(right=0.99, top=0.98, left=0.15, bottom=0.12)
    plt.xlabel(eixo_x_label, fontsize=20, labelpad=8)
    plt.ylabel(eixo_y_label, fontsize=20, labelpad=6)

    ax.tick_params(axis='x', labelsize=16) # Increase the font size of the axis tick labels
    ax.tick_params(axis='y', labelsize=16)

    plt.plot(x1, y1, '#00FF00', label=r'Antiviral $\rho$', linewidth=3)  # Plots the three curves
    plt.plot(x2, y2, '#FFA500', label=r'Antiviral $\xi$', linewidth=3)
    plt.plot(x3, y3, '#0000FF', label=r'Antiviral $\psi$', linewidth=3)

    plt.grid(color='gray', linestyle='--', linewidth=0.3) # Adds grid and legend
    plt.legend(fontsize=18)
    plt.tight_layout()

    final_name = nome + '.png'      # Saves the graph
    save_file(final_name,folder_name)           # Saves the graph using the save_file function
    plt.close(fig)  # Closes the figure to free up memory
    return 0

def graph_generator_from_4_csv(day, csv1, csv2, csv3, csv4, eixo_y_label, nome,eixo_x,eixo_y,folder_name):
    
    data1 = pd.read_csv(csv1)  # Reads data from CSV files
    data2 = pd.read_csv(csv2)
    data3 = pd.read_csv(csv3)
    data4 = pd.read_csv(csv4)

    x1, y1 = data1[eixo_x], data1[eixo_y]  # Assumes the CSV files have 'x' and 'y' columns
    x2, y2 = data2[eixo_x], data2[eixo_y]
    x3, y3 = data3[eixo_x], data3[eixo_y]
    x4, y4 = data4[eixo_x], data4[eixo_y]

    fig, ax = plt.subplots(figsize=(8, 6))  # Creates the figure and axes
    ax.set_facecolor("#FFFFFF")
    eixo_x_label = "Tempo (em dias)"
    plt.subplots_adjust(right=0.99, top=0.98, left=0.15, bottom=0.12)
    plt.xlabel(eixo_x_label, fontsize=20, labelpad=8)
    plt.ylabel(eixo_y_label, fontsize=20, labelpad=1)

    ax.tick_params(axis='x', labelsize=16) # Increase the font size of the axis tick labels
    ax.tick_params(axis='y', labelsize=16)

    plt.plot(x2, y2, '#00FF00', label=r'Antiviral $\rho$', linewidth=3)  # Plots the three curves
    plt.plot(x3, y3, '#FFA500', label=r'Antiviral $\xi$', linewidth=3)
    plt.plot(x4, y4, '#0000FF', label=r'Antiviral $\psi$', linewidth=3)
    plt.plot(x1, y1, "#CC1100", label='Original', linewidth=3)

    plt.axvline(day, color = '#FF00FF', linestyle = '--', linewidth = 1.3,label = 'Início')

    plt.grid(color='gray', linestyle='--', linewidth=0.3) # Adds grid and legend
    plt.legend(fontsize=16)
    plt.tight_layout()

    final_name = nome + '.png'      # Saves the graph
    save_file(final_name,folder_name)           # Saves the graph using the save_file function
    plt.close(fig)  # Closes the figure to free up memory
    return 0

def antiviral_behavior_graph(response, percent_doses, data, label,folder_name):
    fig, ax = plt.subplots(figsize=(8, 6))  # Creates the figure and axes
    ax.set_facecolor("#FFFFFF")
    
    plt.plot(response*100, percent_doses,label=f"{label}",color=data['color'],  linestyle=data["ls"], linewidth=3)  # Plot the dose-response curve
    plt.ylabel("Concentração relativa (%)", fontsize=20)  # Set y-axis label
    plt.xlabel("Eficácia antiviral (em %)", fontsize=20)  # Set x-axis label
    # plt.title("Curva Dose-Resposta: Eficácia antiviral normalizada", fontsize=14)  # (Optional) Set the plot title
    plt.axvline(50, color='#FFD700', linestyle='--', label=f"Eficácia 50% (EC$_{{50}}$ ≈ {data['ec50']*1000:.1f} µg/mL)", linewidth=2)  # Draw a vertical line at 50% efficacy
    
    ax.tick_params(axis='x', labelsize=16) # Increase the font size of the axis tick labels
    ax.tick_params(axis='y', labelsize=16)

    plt.legend(fontsize=16)  # Show the legend
    plt.grid(color='gray', linestyle='--', linewidth=0.3) # Adds grid and legend    
    
    plt.xlim(0, 100)  # Set x-axis limits
    plt.ylim(0, 100)  # Set y-axis limits
    plt.tight_layout()  # Adjust layout to prevent overlap
    
    final_name = "antiviral_behavior_graph.png"  # Name of the output file
    save_file(final_name, folder_name)  # Save the figure using the save_file function
    plt.close(fig)  # Close the figure to free up memory
    return 0 