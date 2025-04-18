# import numpy as np
import pandas as pd
import random as rd
import time
time.time()

from generateGraphs import *
from config_variables import *

for w in range(num_levels):            # Loop through each remediocount value
    for n in range(time - 1):          # Loop through each time step (except the last one)
        for r in range(kutta):         # Loop through the Runge-Kutta stages

            ds[w,n,r] = wks_mu[0] - alpha * s[w,n,r] - a[2] * s[w,n,r] * v[w,n,r]                             # Rate of change for susceptible monocytes
            di[w,n,r] = a[2] * s[w,n,r] * v[w,n,r] - beta * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]               # Rate of change for infected monocytes
            dv[w,n,r] = (1 - remediocount[w]) * k * i[w,n,r] - gamma * v[w,n,r] - a[2] * s[w,n,r] * v[w,n,r]  # Rate of change for viral particles
            dz[w,n,r] = wks_eta[0] + wks_c[0] * i[w,n,r] + wks_d[0] * i[w,n,r] * z[w,n,r] - delta * z[w,n,r]  # Rate of change for T lymphocytes

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

remediocount = [i * 100 for i in remediocount]  # Scale the remediocount values by multiplying each by 100

graph_2D_generator(t, remediocount, s[:, 128, 0], u'Monócitos susceptíveis/\u03bcL', 'gráfico-s')  # Generate a graph for susceptible monocytes (s) over remediocount
graph_2D_generator(t, remediocount, i[:, 128, 0], u'Monócitos infectados/\u03bcL', 'gráfico-i')    # Generate a graph for infected monocytes (i) over remediocount
graph_2D_generator(t, remediocount, v[:, 128, 0], u'Partículas virais/\u03bcL', 'gráfico-v')       # Generate a graph for viral particles (v) over remediocount
graph_2D_generator(t, remediocount, z[:, 128, 0], u'Linfócitos T/\u03bcL', 'gráfico-z')            # Generate a graph for T lymphocytes (z) over remediocount