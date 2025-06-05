from project_libraries import *

## Parameters for the model
a = 0.003        # Successful invasion rate into a susceptible monocyte
mu = 80          # Monocytes produced per day per microliter
eta = 0.265      # T lymphocytes produced per day per microliter for equilibrium of 2000 T cells in absence of infection
alpha = 1/3      # 1/lifespan of a monocyte in days
beta = 1/0.5     # 1/infection period of a monocyte
gamma = 0.8      # Virus release rate
k = 20           # Virus multiplication rate
nu = 0.001       # Elimination rate of infected monocytes
delta = 1/365    # 1/lifespan of T lymphocytes
c = 0.01         # Stimulation of T lymphocyte production by infected monocyte density
d = 0.03         # Stimulation of T lymphocyte production by contact with infected monocytes
beta1 = beta + eta * nu / delta   # Adjusted infection rate
c1 = c + d * eta / delta          # Adjusted stimulation rate

##################

## Antiviral behavior equation
# Logistic (Hill) function                                  # Defines a sigmoid function for dose-response
def sigmoid(x, bottom, top, ec50, hill_slope):              # x: concentration, bottom/top: min/max response, ec50: half-max, hill_slope: steepness
    return bottom + (top - bottom) / (1 + (ec50 / x)**hill_slope)

# Corrected parameters (µg/µL)                              # Stores EC50, Hill slope, color, and linestyle for the compound
params = {
    "Pecto + Acacetina-7-O-rutinosídeo": {
        "ec50": 11.1 / 1000,                               # EC50 value in µg/µL
        "hill": 1.0,                                       # Hill slope
        "color": "#7B68EE",                                # Plot color
        "ls": "-"                                          # Line style
    }
}

dose_max = 0.2                                             # Maximum dose in µg/µL (covers the full curve)

# Relative concentration from 0% to 100%                   # Generates 1001 points from 0.001% to 100% to avoid division by zero
percent_doses = np.linspace(0.001, 100, 1001)
concentrations = (percent_doses / 100) * dose_max          # Converts percent to µg/µL

data = params["Pecto + Acacetina-7-O-rutinosídeo"]         # Retrieves parameter set for the compound
label = "Pecto + Acacetina-7-O-rutinosídeo"                # Label for plotting or reference
response = sigmoid(concentrations, bottom=0, top=1,        # Calculates response curve using the sigmoid function
                   ec50=data["ec50"], hill_slope=data["hill"])


## Other parameters for the antiviral effect (optional)
remediocount_linear = [(i * 0.001) for i in range(1001)]    # Array of linears antiviral levels (from 0 to 1 in steps of 0.001)
remediocount_q = [(i * 0.001)**2 for i in range(1001)]      # Array of quadratic antiviral levels

remediocount_exp = [np.exp(-i * 6* 0.001)**-1 for i in range(1001)] # Array ofexponencial antiviral levels
max_exp = max(remediocount_exp)
remediocount_exp = [val / max_exp for val in remediocount_exp]  # Normalizing to start at 1 and end at 0

## Define the antiviral being used for the model
remedycount = response
num_levels = len(remedycount)                       # Size of the vector that stores the antiviral values

###################################

## Parameters for the Runge-Kutta method
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