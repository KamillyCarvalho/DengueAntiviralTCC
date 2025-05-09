from project_libraries import *  # Importing necessary libraries for numerical calculations and data handling
from generate_graphs import * 
from config_variables import *
# time.time()

antiviral = "rho"  #Options: original; rho; xi; psi  # the original just has to be simulate once   
day = 3           #Options: 1/1.28/2/3/4/5

for w in range(1):
    for n in range(time_counts - 1):          # Loop through each time step (except the last one)
        if(n <= day*100): #100/133/200/300/400/500
            remedy_dose = 0
        else:
            remedy_dose = 1
        for r in range(kutta):                # Loop through the Runge-Kutta stages
            if antiviral == "original":
                ds[w,n,r] = wks_mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
                di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                dv[w,n,r] = k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
                dz[w,n,r] = wks_eta + wks_c * i[w,n,r] + wks_d * i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
            elif antiviral == "rho":
                ds[w,n,r] = wks_mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
                di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                dv[w,n,r] = (1 - remedy_dose) * k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
                dz[w,n,r] = wks_eta + wks_c * i[w,n,r] + wks_d * i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
            elif antiviral == "xi":
                ds[w,n,r] = wks_mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r] * (1 - remedy_dose)
                di[w,n,r] = (1 - remedy_dose) * a * s[w,n,r] * v[w,n,r] - beta * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                dv[w,n,r] = k * i[w,n,r] - gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r] * (1 - remedy_dose)
                dz[w,n,r] = wks_eta + wks_c * i[w,n,r] + wks_d * i[w,n,r] * z[w,n,r] - delta * z[w,n,r]
            elif antiviral == "psi":
                ds[w,n,r] = wks_mu - alpha * s[w,n,r] - a * s[w,n,r] * v[w,n,r]
                di[w,n,r] = a * s[w,n,r] * v[w,n,r] - beta * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
                dv[w,n,r] = k * i[w,n,r] - (1 + remedy_dose) * gamma * v[w,n,r] - a * s[w,n,r] * v[w,n,r]
                dz[w,n,r] = wks_eta + wks_c * i[w,n,r] + wks_d * i[w,n,r] * z[w,n,r] - delta * z[w,n,r]

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

save_data_to_csv(t, s[0,:,0], "antiviral_"+antiviral+"_time_s_day"+day, "Tempo", "Monocitos susceptiveis/uL")
save_data_to_csv(t, i[0,:,0], "antiviral_"+antiviral+"_time_i_day"+day, "Tempo", "Monocitos infectados/uL")
save_data_to_csv(t, v[0,:,0], "antiviral_"+antiviral+"_time_v_day"+day, "Tempo", "Particulas virais/uL")
save_data_to_csv(t, z[0,:,0], "antiviral_"+antiviral+"_time_z_day"+day, "Tempo", "Linfocitos T/uL")