from project_libraries import *  # Importing necessary libraries for numerical calculations and data handling
from config_variables import day

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

def graph_generator_3(eixo_x,eixo_y,t,dados,dados2,nome,folder_name):
    day = 1 #modificar 
    fig, ax = plt.subplots()                         # Creates a figure and an axis
    ax.set_facecolor("#FFFFFF")                      # Sets the axis background to white
    plt.subplots_adjust(right=0.97, top=0.97)        # Adjusts the position of the subplots
    plt.xlabel(eixo_x, fontsize = 12, labelpad = 4)  # Sets the x-axis label
    plt.ylabel(eixo_y,  fontsize = 12, labelpad = 6) # Sets the y-axis label
    plt.plot(t,dados[0,:,0],'#00c2b0',label = 'Com antiviral',linewidth = '1.9')  # Plots the data with antiviral
    plt.plot(t,dados2[0,:,0],'#ff8201',label = 'Sem antiviral',linewidth = '1.9') # Plots the data without antiviral
    plt.axvline(day, color = 'b', linestyle = '--', linewidth = 1.3,label = 'Adm. antiviral') # Adds a vertical line indicating the day of antiviral administration
    leg = ax.legend()                                # Creates the legend
    plt.grid(color = 'gray', linestyle = '-', linewidth = 0.15) # Adds a grid to the graph
    final_name = nome + '[' + str(day) + '].png'     # Final file name, including the day of antiviral administration
    save_file(final_name,folder_name)                            # Saves the graph using the save_file function
    return 0

def graph_2D_generator(t,data1, data2, eixo_y_label, nome,folder_name):
    fig, ax = plt.subplots()                               # Creates a figure and an axis
    ax.set_facecolor("#FFFFFF")                            # Sets the axis background to white
    eixo_x_label = "Antiviral (em %)"                      # Sets the x-axis label
    plt.subplots_adjust(right=0.97, top=0.97)              # Adjusts the position of the subplots
    plt.xlabel(eixo_x_label, fontsize=12, labelpad=4)      # Sets the x-axis label
    plt.ylabel(eixo_y_label, fontsize=12, labelpad=6)      # Sets the y-axis label
    plt.plot(data1, data2, '#00c2b0', label='Dados', linewidth=1.9)  # Plots the data
    plt.grid(color='gray', linestyle='-', linewidth=0.15)  # Adds a grid to the graph
    final_name = nome + '.png'                             # Final file name
    save_file(final_name,folder_name)                                  # Saves the graph using the save_file function
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
    eixo_x_label = "Antiviral (em %)"
    plt.subplots_adjust(right=0.99, top=0.98)
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
    save_file(final_name,folder_name)           # Saves the graph using the save_file function
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

    # Aumenta o tamanho da fonte dos valores nos eixos
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    plt.plot(x1, y1, '#FF0000', label='Original', linewidth=3)
    plt.plot(x2, y2, '#00FF00', label='Antiviral 1', linewidth=3)  # Plots the three curves
    plt.plot(x3, y3, '#FFA500', label='Antiviral 2', linewidth=3)
    plt.plot(x4, y4, '#0000FF', label='Antiviral 3', linewidth=3)
    
    plt.axvline(day, color = '#FF00FF', linestyle = '--', linewidth = 1.3,label = 'Início')

    plt.grid(color='gray', linestyle='--', linewidth=0.3) # Adds grid and legend
    plt.legend(fontsize=16)

    final_name = nome + '.png'      # Saves the graph
    save_file(final_name,folder_name)           # Saves the graph using the save_file function
    return 0