from project_libraries import *  # Importing necessary libraries for numerical calculations and data handling


def save_file(nome):                       # Saves the file in a specific results folder
    output_dir = "resultados"              # Name of the output folder
    if not os.path.exists(output_dir):     # Checks if the folder exists
        os.makedirs(output_dir)            # Creates the folder if it does not exist
    name = os.path.join(output_dir, nome)  # Full path of the file
    plt.savefig(name, format='png')        # Saves the graph in PNG format
    return 0

def graph_generator_3(eixo_x,eixo_y,t,dados,dados2,nome):
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
    save_file(final_name)                            # Saves the graph using the save_file function
    return 0

def graph_2D_generator(t,data1, data2, eixo_y_label, nome):
    fig, ax = plt.subplots()                               # Creates a figure and an axis
    ax.set_facecolor("#FFFFFF")                            # Sets the axis background to white
    eixo_x_label = "Antiviral (em %)"                      # Sets the x-axis label
    plt.subplots_adjust(right=0.97, top=0.97)              # Adjusts the position of the subplots
    plt.xlabel(eixo_x_label, fontsize=12, labelpad=4)      # Sets the x-axis label
    plt.ylabel(eixo_y_label, fontsize=12, labelpad=6)      # Sets the y-axis label
    plt.plot(data1, data2, '#00c2b0', label='Dados', linewidth=1.9)  # Plots the data
    plt.grid(color='gray', linestyle='-', linewidth=0.15)  # Adds a grid to the graph
    final_name = nome + '.png'                             # Final file name
    save_file(final_name)                                  # Saves the graph using the save_file function
    return 0