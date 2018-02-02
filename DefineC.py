#########################################################################
#                                                                       #
#  Create the file c.txt for the BEM.py code                            #
#                                                                       #
#  Alija Bajramovic                                                     #
#  Valentin Bernard                                                     #
#  Naven Goutham                                                        #
#                                                                       #
#########################################################################

# Lbraries #############################################################################

import math
import numpy as np
import matplotlib.pyplot as plt

# Functions ############################################################################

def LiftCoeff (alpha):
    if alpha<21:
        Cl = 0.42625 + 0.11628 * alpha - 0.00063973 * alpha**2 - 8.712 * 10**(-5) * alpha**3 - 4.2576 * 10**(-6)*alpha**4
    else:
        Cl = 0.95
    return Cl

# Definition of Variales ###############################################################

# Physical Parameters
U0 = 10     # Incoming wind speed in m/s
R = 50      # length of the blade in m
N = 3       # Number of blades
Lambda = 8  # Tip speed ratio
Omega = Lambda*U0/R     # Rotor angular velocity in radiants/s
Cl = LiftCoeff(14)

# Computational Parmeters
n = 100       # Number of blade elements
dr = R/100    # Blade element length
r = np.arange(dr, R+dr, dr) - dr/2  # Mid point radius of every element

# Results
Phi = np.zeros(n)   # Angle of incoming wind in radiants
c = np.zeros(n)     # Chord at each element

# Main Program #######################################################################

file = open("c.txt","w") 
for i in range(n):
    # Calculationg the angle of the incoming wind
    Phi[i] = math.atan(U0*(2/3)/(Omega*r[i]))
    c[i] = 8*math.pi*r[i]/N/Cl*(1-math.cos(Phi[i]))
    file.write(str(c[i]) + "\n")
file.close() 

# End of the iteration

plt.figure()
plt.plot(r, c, label = "c")
plt.legend()
plt.show()
