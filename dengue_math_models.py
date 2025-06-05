from project_libraries import *  # Importing necessary libraries for numerical calculations and data handling
from generate_graphs import * 
from config_variables import *

# time.time()

def model_1(antiviral,folder_name_for_figures,folder_name_for_csv):
    for w in range(num_levels):                   # Loop through each remediocount value
        for n in range(time_counts - 1):          # Loop through each time step (except the last one)
            for r in range(kutta):                # Loop through the Runge-Kutta stages
                if antiviral == "rho":
                    ds[w,n,r] = mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
                    di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta1 * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                    dv[w,n,r] = (1 - remedycount[w]) * k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
                    dz[w,n,r] = eta + c1* i[w,n,r] + d* i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
                elif antiviral == "xi":
                    ds[w,n,r] = mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r] * (1 - remedycount[w])
                    di[w,n,r] = (1 - remedycount[w]) * a * s[w,n,r] * v[w,n,r] - beta1 * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                    dv[w,n,r] = k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r] * (1 - remedycount[w])
                    dz[w,n,r] = eta + c1* i[w,n,r] + d* i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
                elif antiviral == "psi":
                    ds[w,n,r] = mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
                    di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta1 * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                    dv[w,n,r] = k * i[w,n,r] - (1 + remedycount[w]) * gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
                    dz[w,n,r] = eta + c1* i[w,n,r] + d* i[w,n,r] * z[w,n,r] - delta * z[w,n,r]

                if r < (kutta - 1):  # Update intermediate stages for Runge-Kutta
                    s[w,n,r+1] = s[w,n,r] + h1 * ds[w,n,r]
                    i[w,n,r+1] = i[w,n,r] + h1 * di[w,n,r]
                    v[w,n,r+1] = v[w,n,r] + h1 * dv[w,n,r]
                    z[w,n,r+1] = z[w,n,r] + h1 * dz[w,n,r]

            # Increment for S, I, V e T
            incremento_s[w,n] = (h/6) * (ds[w,n,0] + 2 * ds[w,n,1] + 2 * ds[w,n,2] + ds[w,n,3])  
            incremento_i[w,n] = (h/6) * (di[w,n,0] + 2 * di[w,n,1] + 2 * di[w,n,2] + di[w,n,3])  
            incremento_v[w,n] = (h/6) * (dv[w,n,0] + 2 * dv[w,n,1] + 2 * dv[w,n,2] + dv[w,n,3])  
            incremento_z[w,n] = (h/6) * (dz[w,n,0] + 2 * dz[w,n,1] + 2 * dz[w,n,2] + dz[w,n,3])  
                
            # Update S, I, V e T for the next time step
            s[w,n+1,0] = s[w,n,0] + incremento_s[w,n]  
            i[w,n+1,0] = i[w,n,0] + incremento_i[w,n]
            v[w,n+1,0] = v[w,n,0] + incremento_v[w,n]
            z[w,n+1,0] = z[w,n,0] + incremento_z[w,n]

    remediocount = [i * 100 for i in remedycount]  # Scale the remediocount values by multiplying each by 100

    # use to debug 
    # graph_2D_generator(t, remediocount, s[:, 128, 0], u'Monócitos susceptíveis/\u03bcL', 'grafico-s',folder_name_for_figures)  # Generate a graph for susceptible monocytes (s) over remediocount
    # graph_2D_generator(t, remediocount, i[:, 128, 0], u'Monócitos infectados/\u03bcL', 'grafico-i',folder_name_for_figures)    # Generate a graph for infected monocytes (i) over remediocount
    # graph_2D_generator(t, remediocount, v[:, 128, 0], u'Partículas virais/\u03bcL', 'grafico-v',folder_name_for_figures)       # Generate a graph for viral particles (v) over remediocount
    # graph_2D_generator(t, remediocount, z[:, 128, 0], u'Linfócitos T/\u03bcL', 'grafico-z',folder_name_for_figures)            # Generate a graph for T lymphocytes (z) over remediocount

    save_data_to_csv(remediocount, s[:, 127, 0], "susceptible_monocytes_antiviral_"+antiviral, "Antiviral", "Monocitos susceptiveis/uL",folder_name_for_csv)
    save_data_to_csv(remediocount, i[:, 127, 0], "infected_monocytes_antiviral_"+antiviral, "Antiviral", "Monocitos infectados/uL",folder_name_for_csv)
    save_data_to_csv(remediocount, v[:, 127, 0], "viral_particles_antiviral_"+antiviral, "Antiviral", "Particulas virais/uL",folder_name_for_csv)
    save_data_to_csv(remediocount, z[:, 127, 0], "t_lymphocytes_antiviral_"+antiviral, "Antiviral", "Linfocitos T/uL",folder_name_for_csv)

def model_2(antiviral,folder_name_for_figures,folder_name_for_csv,day):
    w=0
    for n in range(time_counts - 1):          # Loop through each time step (except the last one)
        if(n <= day*100): #100/128/200/300/400/500
            remedy_dose = 0
        else:
            remedy_dose = remedycount[int(n-(day*100))]
        for r in range(kutta):                # Loop through the Runge-Kutta stages
            if antiviral == "rho":
                ds[w,n,r] = mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
                di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta1 * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                dv[w,n,r] = (1 - remedy_dose) * k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
                dz[w,n,r] = eta + c1* i[w,n,r] + d* i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
            elif antiviral == "xi":
                ds[w,n,r] = mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r] * (1 - remedy_dose)
                di[w,n,r] = (1 - remedy_dose) * a * s[w,n,r] * v[w,n,r] - beta1 * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                dv[w,n,r] = k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r] * (1 - remedy_dose)
                dz[w,n,r] = eta + c1* i[w,n,r] + d* i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
            elif antiviral == "psi":
                ds[w,n,r] = mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
                di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta1 * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                dv[w,n,r] = k * i[w,n,r] - (1 + remedy_dose) * gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
                dz[w,n,r] = eta + c1* i[w,n,r] + d* i[w,n,r] * z[w,n,r] - delta * z[w,n,r]

            if r < (kutta - 1):  # Update intermediate stages for Runge-Kutta
                s[w,n,r+1] = s[w,n,r] + h1 * ds[w,n,r]
                i[w,n,r+1] = i[w,n,r] + h1 * di[w,n,r]
                v[w,n,r+1] = v[w,n,r] + h1 * dv[w,n,r]
                z[w,n,r+1] = z[w,n,r] + h1 * dz[w,n,r]

        # Increment for S, I, V e T
        incremento_s[w,n] = (h/6) * (ds[w,n,0] + 2 * ds[w,n,1] + 2 * ds[w,n,2] + ds[w,n,3])  
        incremento_i[w,n] = (h/6) * (di[w,n,0] + 2 * di[w,n,1] + 2 * di[w,n,2] + di[w,n,3])  
        incremento_v[w,n] = (h/6) * (dv[w,n,0] + 2 * dv[w,n,1] + 2 * dv[w,n,2] + dv[w,n,3])  
        incremento_z[w,n] = (h/6) * (dz[w,n,0] + 2 * dz[w,n,1] + 2 * dz[w,n,2] + dz[w,n,3])  
        
        # Update S, I, V e T for the next time step
        s[w,n+1,0] = s[w,n,0] + incremento_s[w,n]  
        i[w,n+1,0] = i[w,n,0] + incremento_i[w,n]
        v[w,n+1,0] = v[w,n,0] + incremento_v[w,n]
        z[w,n+1,0] = z[w,n,0] + incremento_z[w,n]

    save_data_to_csv(t, s[0,:,0], "antiviral_"+antiviral+"_time_s_day"+str(day), "Tempo", "Monocitos susceptiveis/uL",folder_name_for_csv)
    save_data_to_csv(t, i[0,:,0], "antiviral_"+antiviral+"_time_i_day"+str(day), "Tempo", "Monocitos infectados/uL",folder_name_for_csv)
    save_data_to_csv(t, v[0,:,0], "antiviral_"+antiviral+"_time_v_day"+str(day), "Tempo", "Particulas virais/uL",folder_name_for_csv)
    save_data_to_csv(t, z[0,:,0], "antiviral_"+antiviral+"_time_z_day"+str(day), "Tempo", "Linfocitos T/uL",folder_name_for_csv)

def model_3(antiviral,folder_name_for_figures,folder_name_for_csv):
    w=0
    for n in range(time_counts - 1):          # Loop through each time step (except the last one)
        for r in range(kutta):                # Loop through the Runge-Kutta stages
            ds[w,n,r] = mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
            di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta1 * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
            dv[w,n,r] = k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
            dz[w,n,r] = eta + c1* i[w,n,r] + d* i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
            
            if r < (kutta - 1):  # Update intermediate stages for Runge-Kutta
                s[w,n,r+1] = s[w,n,r] + h1 * ds[w,n,r]
                i[w,n,r+1] = i[w,n,r] + h1 * di[w,n,r]
                v[w,n,r+1] = v[w,n,r] + h1 * dv[w,n,r]
                z[w,n,r+1] = z[w,n,r] + h1 * dz[w,n,r]

        # Increment for S, I, V e T
        incremento_s[w,n] = (h/6) * (ds[w,n,0] + 2 * ds[w,n,1] + 2 * ds[w,n,2] + ds[w,n,3])  
        incremento_i[w,n] = (h/6) * (di[w,n,0] + 2 * di[w,n,1] + 2 * di[w,n,2] + di[w,n,3])  
        incremento_v[w,n] = (h/6) * (dv[w,n,0] + 2 * dv[w,n,1] + 2 * dv[w,n,2] + dv[w,n,3])  
        incremento_z[w,n] = (h/6) * (dz[w,n,0] + 2 * dz[w,n,1] + 2 * dz[w,n,2] + dz[w,n,3])  
        
        # Update S, I, V e T for the next time step
        s[w,n+1,0] = s[w,n,0] + incremento_s[w,n]  
        i[w,n+1,0] = i[w,n,0] + incremento_i[w,n]
        v[w,n+1,0] = v[w,n,0] + incremento_v[w,n]
        z[w,n+1,0] = z[w,n,0] + incremento_z[w,n]

    save_data_to_csv(t, s[0,:,0], "antiviral_"+antiviral+"_time_s", "Tempo", "Monocitos susceptiveis/uL",folder_name_for_csv)
    save_data_to_csv(t, i[0,:,0], "antiviral_"+antiviral+"_time_i", "Tempo", "Monocitos infectados/uL",folder_name_for_csv)
    save_data_to_csv(t, v[0,:,0], "antiviral_"+antiviral+"_time_v", "Tempo", "Particulas virais/uL",folder_name_for_csv)
    save_data_to_csv(t, z[0,:,0], "antiviral_"+antiviral+"_time_z", "Tempo", "Linfocitos T/uL",folder_name_for_csv)

    graph_2D_generator(t, t, s[0, :, 0], u'Monócitos susceptíveis/\u03bcL', 'graph-s-original',folder_name_for_figures)  # Generate a graph for susceptible monocytes (s) over remediocount
    graph_2D_generator(t, t, i[0, :, 0], u'Monócitos infectados/\u03bcL', 'graph-i-original',folder_name_for_figures)    # Generate a graph for infected monocytes (i) over remediocount
    graph_2D_generator(t, t, v[0, :, 0], u'Partículas virais/\u03bcL', 'graph-v-original',folder_name_for_figures)       # Generate a graph for viral particles (v) over remediocount
    graph_2D_generator(t, t, z[0, :, 0], u'Linfócitos T/\u03bcL', 'graph-z-original',folder_name_for_figures)            # Generate a graph for T lymphocytes (z) over remediocount

    # print(np.argmax(v[0,:,0]))