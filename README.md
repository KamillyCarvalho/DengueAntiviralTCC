# DengueAntiviralTCC

This project aims to simulate and analyze the effects of antiviral treatments on dengue infection dynamics. Using mathematical modeling and numerical methods, the project generates data and visualizations to better understand the behavior of health and infected monocytes, viral particles, and T lymphocytes under different antiviral scenarios.

The repository is organized into modular Python scripts for better maintainability and scalability. Each script has a specific purpose, such as managing dependencies, defining configuration variables, generating graphs, and processing data.

### Key Features:
- Simulation of dengue infection dynamics using Runge-Kutta numerical methods.
- Visualization of results through customizable graphs.
- Modular structure for easy maintenance and extension.
- Data export functionality for further analysis.

## How to Run the `main` Function
antiviral = "rho"  # Options: "rho", "xi", "psi"

## project_libraries.py

The `project_libraries.py` script is responsible for managing and organizing the libraries or dependencies used in this project. 

### Key Features:
- Listing all required libraries for the project.
- Checking if the necessary libraries are installed.
- Avoiding the repetition of library imports across multiple files by centralizing them in one location.

## config_variables.py

The `config_variables.py` file is responsible for managing configuration settings used throughout the project. It typically contains constants, environment variables, or other settings that can be easily modified to adapt the application to different environments (e.g., development, testing, production).

### Key Features:
- Centralized location for configuration values.
- Simplifies the process of updating settings without modifying the main application code.

## generate_graphs.py

The `generateGraphs.py` file contains functions to generate and save graphs for the Dengue Antiviral Project.

### Functions

- **`save_file(nome)`**  
    Saves the graph as a PNG file in the `results` folder.

- **`save_data_to_csv(x_data, y_data, file_name, x_label, y_label)`**  
    Saves the data to a csv file in the `results` folder.

- **`graph_generator_3(eixo_x, eixo_y, t, dados, dados2, nome)`**  
    Creates a graph comparing two datasets with a vertical line indicating the day of antiviral administration.

- **`graph_2D_generator(t, data1, data2, eixo_y_label, nome)`**  
    Generates a 2D graph with a fixed x-axis label ("Antiviral (em %)").