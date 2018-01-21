#########################################################################
#                                                                       #
#  Performance analysis of a wind turbine for a given incoming wind     #
#                                                                       #
#  Version 1, 16/01/2018                                                #
#  Valentin Bernard                                                     #
#                                                                       #
#########################################################################

# Lbraries #############################################################################

import math
import matplotlib.pyplot as plt
import numpy as np


# Functions ############################################################################

def finiteSums( function, dx ):
    "Summs the values in an array function multiplied by a stepsize dx"
    integral = 0
    for value in function:
        integral += value*dx
    return value

# Definition of Variales ###############################################################

# Physical Parameters
U0 = 10       # Incoming wind speed in m/s
Rho = 1       # Air density in kg/m^3
R = 50        # length of the blade in m
N = 3         # Number of blades
Lambda = 13   # Tip speed ratio

# Computational Parmeters
n = 100            # Number of blade elements
dx = R/100         # Blade element length

# Blade characteristics

cd = np.genfromtxt("LiftCoeff.txt", unpack=True)  # Drag coefficient as f(alpha), 1 value per degree
cl = [0] * 30       # Drag coefficient as f(alpha), 1 value per degree
c = [0] * n        # Airfoil chord in m as function of r
Beta = [0] * n     # Twist angle in radiants as function of r

# Derived quantities
Omega = Lambda*U0/R     # Rotor angular velocity in radiants/s
Area = R*R*math.pi      # Rotor area in m^2
Pin = math.pow(U0,3)*Rho*Area     # Incoming wind power in W


# Results
Cp = 0             # Power Coefficient
Ct = 0             # Thrust coefficient
dL = [0] * n       # Lift on each elment in N
dD = [0] * n       # Drag on each elment in N
dM = [0] * n       # Momentm on each elment in Nm
M = 0              # Total momntum in Nm
P = 0              # Power in W
dT = [0] * n       # Thrust force on each element in N
T = 0              # Thrust force
Phi = [0] * n      # Angle of incoming wind in radiants
Alpha = [0] * n    # Angle of attack in radiants

# Main Program #######################################################################

a = 0              # Linear induction factor
aa = 0             # Angular induction factor
U = U0*(1-a)       # Wind speed at the turbine

for i in range(n):
    Phi[i] = math.atan(U/Omega/dx/(i+1)/(1-aa))*180/math.pi
    Alpha[i] = Phi[i]-Beta[i]
    dL[i] = cl[i]/2*Rho*U*U*c[i]*dx
    dD[i] = cd[i]/2*Rho*U*U*c[i]*dx
    
    
    

M = sum(dM)        # Calculation of total momentum
P = M*Omega        # Calculation of output power
Cp = P/Pin         # Calculation of Power coefficient


