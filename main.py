# import numpy as np
import pandas as pd
import random as rd
import time
from mpl_toolkits.mplot3d import Axes3D
time.time()

from generateGraphs import *

wks_s = np.array([1.0, 1.96, 2.425, 3.63])           #% da população de monócitos ao longo da gravidez em relação a não grávida
wks_z = np.array([1.0, 0.81, 0.64, 0.57])            #% da população de linfócitos ao longo da gravidez em relação a não grávida
wks_IFN = np.array([1.0, 0.5963, 0.5611, 1.0704])    #% da produção de linfócitos-T com o mesmo estímulo
a = np.array([0.001, 0.002, 0.003])                  #taxa de invasão bem-sucedida em um monócito susceptível
wks_mu = 80 * wks_s                                  #monócitos produzidos/dia.uL
wks_eta = 0.265 * wks_z                              #linfo-T produzidos/dia.uL para equilíbrio de 2000 linfo-T na ausência de infecção
alpha = 1/3                                          #1/período de vida de um monócito em dias
beta = 1/0.5                                         #1/período de infecção de um monócito
gamma = 0.8                                          #taxa de liberação de vírus
k = 20                                               #taxa de multiplicação de vírus
nu = 0.001                                           #taxa de eliminação de monócito infectado
delta = 1/365                                        #1/período de vida de linfo-T
wks_c = 0.01 * wks_IFN                               #estímulo de produção de linfócitos-T pela densidade de monócitos infectados
wks_d = 0.03 * wks_IFN                               #estímulo de produção de linfócitos-T pelos contatos com monócitos infectados
wks_beta1 = beta + wks_eta * nu / delta
wks_c1 = wks_c + wks_d * wks_eta / delta

h = 0.01                                             #parâmetro de precisão para método numérico de Runge-Kutta
h1 = h/2                                             #parâmetro intermediário de precisão para método numérico de Runge-Kutta
t = np.arange(0, 10+h, h)                            #vetor de tempo, até 10 dias
time = len(t)
kutta = 4

remediocount = [i * 0.001 for i in range(1001)]
num_levels = len(remediocount)

s = np.zeros((num_levels,time,kutta))       # Array to store susceptible monocytes at each time step and Runge-Kutta stage
v = np.zeros((num_levels,time,kutta))       # Array to store viral particles at each time step and Runge-Kutta stage
i = np.zeros((num_levels,time,kutta))       # Array to store infected monocytes at each time step and Runge-Kutta stage
z = np.zeros((num_levels,time,kutta))       # Array to store T lymphocytes at each time step and Runge-Kutta stage

ds = np.zeros((num_levels,time - 1,kutta))  # Array to store rate of change for susceptible monocytes
dv = np.zeros((num_levels,time - 1,kutta))  # Array to store rate of change for viral particles
di = np.zeros((num_levels,time - 1,kutta))  # Array to store rate of change for infected monocytes
dz = np.zeros((num_levels,time - 1,kutta))  # Array to store rate of change for T lymphocytes

incremento_s = np.zeros((num_levels,time))  # Array to store increments for susceptible monocytes
incremento_i = np.zeros((num_levels,time))  # Array to store increments for infected monocytes
incremento_v = np.zeros((num_levels,time))  # Array to store increments for viral particles
incremento_z = np.zeros((num_levels,time))  # Array to store increments for T lymphocytes

###############################

# Initial values for the populations
s[:,0,0] = 250        # Susceptible monocytes per microliter
i[:,0,0] = 10         # Infected monocytes per microliter
v[:,0,0] = 165        # Viral particles per microliter
z[:,0,0] = 2000       # T lymphocytes per microliter

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