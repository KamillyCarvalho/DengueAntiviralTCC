# DengueAntiviralTCC

## Graph Generation Script

The `generateGraphs.py` file contains functions to generate and save graphs for the Dengue Antiviral Project.

### Functions

- **`Save_file(nome)`**  
    Saves the graph as a PNG file in the `resultados` folder.

- **`graph_generator_3(eixo_x, eixo_y, t, dados, dados2, nome)`**  
    Creates a graph comparing two datasets with a vertical line indicating the day of antiviral administration.

- **`graph_2D_generator(t, data1, data2, eixo_y_label, nome)`**  
    Generates a 2D graph with a fixed x-axis label ("Antiviral (em %)").

### Output

All graphs are saved as PNG files in the `resultados` folder.