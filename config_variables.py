from project_libraries import *

day = 3           #Options: 1/1.28/2/3/4/5

# Parameters for the model
a = 0.003                                            #taxa de invasão bem-sucedida em um monócito susceptível
wks_mu = 80                                          #monócitos produzidos/dia.uL
wks_eta = 0.265                                      #linfo-T produzidos/dia.uL para equilíbrio de 2000 linfo-T na ausência de infecção
alpha = 1/3                                          #1/período de vida de um monócito em dias
beta = 1/0.5                                         #1/período de infecção de um monócito
gamma = 0.8                                          #taxa de liberação de vírus
k = 20                                               #taxa de multiplicação de vírus
nu = 0.001                                           #taxa de eliminação de monócito infectado
delta = 1/365                                        #1/período de vida de linfo-T
wks_c = 0.01                                         #estímulo de produção de linfócitos-T pela densidade de monócitos infectados
wks_d = 0.03                                         #estímulo de produção de linfócitos-T pelos contatos com monócitos infectados
wks_beta1 = beta + wks_eta * nu / delta
wks_c1 = wks_c + wks_d * wks_eta / delta

# Parameters for the antiviral effect
remediocount_linear = [(i * 0.001) for i in range(1001)]    # Array of linears antiviral levels (from 0 to 1 in steps of 0.001)
remediocount_q = [(i * 0.001)**2 for i in range(1001)]      # Array of quadratic antiviral levels

remediocount_exp = [np.exp(-i * 6* 0.001)**-1 for i in range(1001)] # Array ofexponencial antiviral levels
max_exp = max(remediocount_exp)
remediocount_exp = [val / max_exp for val in remediocount_exp]  # Normalizing to start at 1 and end at 0

num_levels = len(remediocount_linear)                       # Size of the vector that stores the antiviral values

# Parameters for the Runge-Kutta method
h = 0.01                                             # Precision parameter for the Runge-Kutta numerical method
h1 = h/2                                             # Intermediate precision parameter for the Runge-Kutta numerical method
t = np.arange(0, 10+h, h)                            # Time vector, to count up to 10 days
time_counts = len(t)                                 # Number of time steps based on the time vector
kutta = 4                                            # Number of stages in the Runge-Kutta method

###############################

# Arrays to store the results of the Runge-Kutta method
s = np.zeros((num_levels,time_counts,kutta))       # Array to store susceptible monocytes at each time step and Runge-Kutta stage
v = np.zeros((num_levels,time_counts,kutta))       # Array to store viral particles at each time step and Runge-Kutta stage
i = np.zeros((num_levels,time_counts,kutta))       # Array to store infected monocytes at each time step and Runge-Kutta stage
z = np.zeros((num_levels,time_counts,kutta))       # Array to store T lymphocytes at each time step and Runge-Kutta stage

ds = np.zeros((num_levels,time_counts - 1,kutta))  # Array to store rate of change for susceptible monocytes
dv = np.zeros((num_levels,time_counts - 1,kutta))  # Array to store rate of change for viral particles
di = np.zeros((num_levels,time_counts - 1,kutta))  # Array to store rate of change for infected monocytes
dz = np.zeros((num_levels,time_counts - 1,kutta))  # Array to store rate of change for T lymphocytes

incremento_s = np.zeros((num_levels,time_counts))  # Array to store increments for susceptible monocytes
incremento_i = np.zeros((num_levels,time_counts))  # Array to store increments for infected monocytes
incremento_v = np.zeros((num_levels,time_counts))  # Array to store increments for viral particles
incremento_z = np.zeros((num_levels,time_counts))  # Array to store increments for T lymphocytes

###############################

# Initial values for the populations
s[:,0,0] = 250        # Susceptible monocytes per microliter
i[:,0,0] = 10         # Infected monocytes per microliter
v[:,0,0] = 165        # Viral particles per microliter
z[:,0,0] = 2000       # T lymphocytes per microliter