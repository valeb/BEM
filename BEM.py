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
import tkinter.filedialog

# Functions ############################################################################


# Definition of Variales ###############################################################

# Physical Parameters
U0 = 10       # Incoming wind speed in m/s
Rho = 1.225       # Air density in kg/m^3
R = 50        # length of the blade in m
N = 3         # Number of blades
Lambda = 7   # Tip speed ratio

# Computational Parmeters
n = 100       # Number of blade elements
dx = R/100    # Blade element length

# Blade characteristics

cl_data = np.genfromtxt("LiftCoeff.txt") # Lift coefficient as f(alpha), 1 value per degree
cd_data = np.genfromtxt("DragCoeff.txt") # Drag coefficient as f(alpha), 1 value per degree
c = [2.0] * n                       # Airfoil chord in m as function of r
Beta = np.genfromtxt("Beta.txt")*math.pi/180    # Twist angle in radiants as function of r

# Derived quantities
Omega = Lambda*U0/R     # Rotor angular velocity in radiants/s
Area = R**2*math.pi      # Rotor area in m^2
Pin = 1/2*U0**3*Rho*Area     # Incoming wind power in W


# Results
Cp = 0.0            # Power Coefficient
Ct = 0.0            # Thrust coefficient
dL = [0.0] * n      # Lift on each elment in N
dD = [0.0] * n      # Drag on each elment in N
dM = [0.0] * n      # Momentm on each elment in Nm
M = 0.0             # Total momntum in Nm
P = 0.0             # Power in W
dT = [0.0] * n      # Thrust force on each element in N
T = 0.0             # Thrust force
Phi = [0.0] * n     # Angle of incoming wind in radiants
Alpha = [0.0] * n   # Angle of attack in radiants

# Main Program #######################################################################

a = [0.0]*n             # Linear induction factor, initial guess
aa = [0.0]*n            # Angular induction factor, initial guess

# Variables for the iteration
Difference = 999        
a_new = [0.0]*n
aa_new = [0.0]*n
counter = 0
#evolution_a_mean = [0]
# Loop for the iteration
while (Difference > 0.0000001) :
    Difference = 0.0
    counter += 1
    for i in range(n):
        # Calculationg the angle of the incoming wind
        Phi[i] = math.atan(U0*(1-a[i])/(Omega*dx*(i+1)*(1+aa[i])))
        Alpha[i] = Phi[i]-Beta[i]
        # Calculationg relative velocity
        Urel  = U0*(1-a[i])/math.sin(Phi[i]) # math.sqrt((U0*(1-a[i]))**2+(Omega*dx*(i+1)*(1+aa[i]))**2)
        # Calculating drag and lift on the element from  interpolated coefficient data
        cl = 0.5*(cl_data[math.floor(Alpha[i])]+cl_data[math.ceil(Alpha[i])])
        cd = 0.5*(cd_data[math.floor(Alpha[i])]+cd_data[math.ceil(Alpha[i])])
        dL[i] = N*cl*1/2*Rho*Urel**2*c[i]*dx
        dD[i] = N*cd*1/2*Rho*Urel**2*c[i]*dx
        # Calculation Momentum and Thrust from the forces
        dM[i] = (i+1)*dx*(dL[i]*math.sin(Phi[i])-dD[i]*math.cos(Phi[i]))
        dT[i] = dL[i]*math.cos(Phi[i])+dD[i]*math.sin(Phi[i])
        # Calculation of the induction factors
        Sigma = N*c[i]/(2*math.pi*dx*(i+1)) # Local solitity
        a_new[i] = 1/(1+(4*math.sin(Phi[i])**2)/(Sigma*cl*math.cos(Phi[i])))
        aa_new[i] = 1/((4*math.cos(Phi[i])/(Sigma*cl))-1)
        Difference += (a_new[i]-a[i])**2+(aa_new[i]-aa[i])**2
        a[i] = a_new[i]
        aa[i] = aa_new[i]
    #evolution_a_mean.append(sum(a)/n)
    print(counter, "Difference= ", Difference)
# End of the iteration
#plt.plot(evolution_a_mean)
plt.plot(Beta, label = "Beta")
plt.plot(Phi, label = "Phi")
plt.plot(Alpha, label = "Alpha")

M = sum(dM)     # Calculation of total momentum
T = sum(dT)     # Calculation of total thrust     
P = M*Omega     # Calculation of output power
Cp = P/Pin      # Calculation of Power coefficient
print("P= ", P , "\nPin= ", Pin, "\nCp= " , Cp)

'''
plt.plot(dM, "g-", label='Momentum')
plt.plot(dT, "r-", label = 'Thrust')
plt.xlabel("Balde Element")
plt.ylabel("Force in N / Momntum in Nm")
plt.title("Momentum and thrust force on each blade eleent") 
'''
plt.legend()
plt.show()



