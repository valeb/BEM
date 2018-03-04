#########################################################################
#                                                                       #
#  Performance analysis of a wind turbine for a given incoming wind     #
#  using blade element momentum theory                                  #
#                                                                       #
#  Alija Bajramovic                                                     #
#  Valentin Bernard                                                     #
#  Naven Goutham                                                        #
#                                                                       #
#########################################################################

# Lbraries #############################################################################

from sys import exit
import math
import matplotlib.pyplot as plt
import numpy as np
import tkinter.filedialog

# Functions ###########################################################################

def Coefficient (alpha, AlphaData, CoeffData): 
    # Interpolate data to obtain Cl or Cd for given alpha in degrees
    i = 0
    if alpha < AlphaData[0]: return CoeffData[0]
    if alpha > AlphaData[-1]: return CoeffData[-1]
    while not (AlphaData[i] <= alpha and AlphaData[i+1] > alpha):
        i += 1
    coeff = (CoeffData[i+1]*(alpha-AlphaData[i]) + CoeffData[i]*(AlphaData[i+1]-alpha))/(AlphaData[i+1]-AlphaData[i]) 
    return coeff


# Variables ###############################################################

# Physical Parameters
U0_data = np.arange(5, 25.1, 0.1)         # Incoming wind speed in m/s
Rho = 1.225                         # Air density in kg/m^3
R = 50                              # length of the blade in m
N = 3                               # Number of blades
Omega = 1.6                         # Rotor angular velocity to obtain Lambda = 8 at 10m/s
Area = R**2*math.pi                 # Rotor area in m^2

# Computational Parmeters
n = 100                             # Number of blade elements
dr = R/n                            # Blade element length
r = np.arange(dr, R+dr, dr)         # Mid point radius of every element
Limit = 10**(-4)                    # Tolerance for iteration

# Import blade characteristics from file
c = np.genfromtxt("./Results/NACA_UC_c.txt")          # Airfoil chord in m as function of r
Beta = np.genfromtxt("./Results/NACA_UC_Beta.txt")    # Twist angle in radians as function of r
f = open("./Results/NACA_UC_Foils.txt","r") 
Foils = f.read().splitlines() 
f.close()

# Main Program #################################################################

Cp = np.zeros(U0_data.size)

for U0 in U0_data:
    Lambda = R*Omega/U0
    Pin = 1/2*U0**3*Rho*Area            # Incoming wind power in W
    Lambda_r = r*Lambda/R               # Local tip speed ratio
    # Results

    dM = [0.0] * n                      # Momentm on each elment in Nm
    Phi = np.zeros(n)                   # Angle of incoming wind in radiants
    Alpha = np.zeros(n)                 # Angle of attack in radiants
    a = [1/3]*n
    aa = [0]*n    

    for i in range(5,n):
        AlphaData = np.genfromtxt(Foils[i], usecols=0)
        ClData = np.genfromtxt(Foils[i], usecols=1)
        CdData = np.genfromtxt(Foils[i], usecols=2)
        Sigma = N*c[i]/(2*math.pi*r[i])
        # Variables for the iteration
        Difference = 10
        counter = 0;
        # Loop for the iteration
        while Difference > Limit and counter < 100 :
        # Calculationg the angle of the incoming wind
            Phi[i] = math.atan((1-a[i])/(1+aa[i])/Lambda_r[i])
            Alpha[i] = Phi[i]*180/math.pi - Beta[i]
        # Calculating drag and lift coefficiet from empiric equation
            Cl = Coefficient(Alpha[i], AlphaData, ClData)
            Cd = Coefficient(Alpha[i], AlphaData, CdData)
        # Calculation of the new induction factors
            a_new = 1/(1+(2*math.sin(Phi[i]))**2/(Sigma*Cl*math.cos(Phi[i])))
            aa_new = 1/((4*math.cos(Phi[i])/(Sigma*Cl))-1)
            Difference = (a_new-a[i])**2+(aa_new-aa[i])**2
            a[i] = a_new
            aa[i] = aa_new
            counter += 1
        # End of the iteration
        # Calculationg relative velocity
        Urel  = U0*(1-a[i])/math.sin(Phi[i])
        # Calculation Momentum and Thrust from the forces
        dM[i] = N*0.5*Rho*Urel**2*(Cl*math.sin(Phi[i])-Cd*math.cos(Phi[i]))*c[i]*r[i]*dr

    P = sum(dM)*Omega
    Cp[int((U0-5)/0.1)] = P/Pin


# Plots ############################################################################

plt.plot(U0_data, Cp)

# File output ####################################################################

# Export resulting Cp values
f = open("./Results/NACA_UC_Cp.txt","w") 
for i in range(Cp.size):
    f.write(str(U0_data[i]) + "\t" + str(Cp[i]) + "\n")
f.close()

plt.show()

