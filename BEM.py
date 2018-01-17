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

# Definition of Variales ###############################################################

# Computational Parmeters
n = 100       # Number of blade elements

# Physical Parameters
U0 = 10       # Incoming wind speed in m/s
Roh = 1       # Air density in kg/m^3
R = 50        # length of the blade in m
N = 3         # Number of blades
Lambda = 13   # Tip speed ratio
cd = [0] * n       # Drag coefficient as function of alpha
c = [0] * n        # Airfoil chord in m as function of r
Beta = [0] * n     # Twist angle in radiants as function of r
a = 0              # Linear induction factor
aa = 0             # Angular induction factor

# Derived quantities
Omega = Lambda*U0/R     # Rotor angular velocity in radiants/s
Area = R*R*math.pi      # Rotor area in m^2
Pin = math.pow(U0,3)*Roh*Area     # Incoming wind power in W

# Results
Cp = 0             # Power Coefficient
Ct = 0             # Thrust coefficient
dM = [0] * n       # Momentm on each elment in Nm
M = 0              # Total momntum in Nm
P = 0              # Power in W
dT = [0] * n       # Thrust force on each element in N
T = 0              # Thrust force
Phi = [0] * n      # Angle of incoming wind in radiants
Alpha = [0] * n    # Angle of attack in radiants

# Main Program #######################################################################

U = U0*(1-a)       # Wind speed at the turbine

M = 0
for value in dM:   # Calculation of total momentm
    M += value*R/100
P = M*Omega        # Calculation of output power
Cp = P/Pin         # Calculation of Power coefficient

