#########################################################################
#                                                                       #
#  Create the file Beta.txt for the BEM.py code                         #
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

# Definition of Variales ###############################################################

# Physical Parameters
U0 = 10                             # Incoming wind speed in m/s
R = 50                              # Length of the blade in m
R_i = 5                             # Inner radius in m
Lambda = 8                          # Tip speed ratio
Omega = Lambda*U0/R                 # Rotor angular velocity in radiants/s

# Computational Parmeters
n = 100                             # Number of blade elements
dr = (R-R_i)/n                      # Blade element length
r = np.arange(dr+R_i, R+dr, dr)     # Mid point radius of every element
Lambda_r = r*Lambda/R               # Local tip speed ratio

# Results
Phi = np.zeros(n)   # Angle of incoming wind in radiants
Beta = np.zeros(n) 

# Main Program #######################################################################

file = open("Beta.txt","w") 
for i in range(n):
    # Calculationg the angle of the incoming wind
    Phi[i] = (2/3)*math.atan(1/Lambda_r[i])
    Beta[i] = Phi[i] - 14/180*math.pi
    file.write(str(Beta[i]) + "\n")
file.close() 

# End of the iteration

plt.figure()
plt.plot(r, Beta/math.pi*180, label = "Beta")
plt.legend()
plt.show()
