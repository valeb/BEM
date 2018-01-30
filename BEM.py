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
Rho = 1225       # Air density in kg/m^3
R = 50        # length of the blade in m
N = 3         # Number of blades
Lambda = 8   # Tip speed ratio

# Computational Parmeters
n = 100       # Number of blade elements
dr = R/100    # Blade element length
r = np.arange(dr, R+dr, dr)     # radius of every element

# Blade characteristics

cl_data = np.genfromtxt("LiftCoeff.txt")        # Lift coefficient as f(alpha), 1 value per degree
cd_data = np.genfromtxt("DragCoeff.txt")        # Drag coefficient as f(alpha), 1 value per degree
c = [1] * n                                   # Airfoil chord in m as function of r
Beta = np.genfromtxt("Beta.txt")*math.pi/180    # Twist angle in radiants as function of r

plt.figure()
plt.plot(cl_data, label = "Lift Coefficient")
plt.plot(cd_data, label = "Drag Coefficient")
plt.xlabel("Angle of Attack in Degrees")
plt.legend()


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
Phi = np.empty(n)   # Angle of incoming wind in radiants
Alpha = np.empty(n) # Angle of attack in radiants

# Main Program #######################################################################

a = [0.0]*n             # Linear induction factor, initial guess
aa = [0.0]*n            # Angular induction factor, initial guess

# Variables for the iteration
Difference = 999        
a_new = [0.0]*n
aa_new = [0.0]*n
counter = 0
evolution_a_mean = [0]
evolution_aa_mean = [0]

# Loop for the iteration
while (Difference > 10**(-3)) :
    Difference = 0.0
    counter += 1
    for i in range(n):
        # Calculationg the angle of the incoming wind
        Phi[i] = math.atan(U0*(1-a[i])/(Omega*dr*(i+1)*(1+aa[i])))
        if Phi[i] > math.pi/2 :
            Phi[i] = Phi[i] - math.pi
        Alpha[i] = Phi[i]-Beta[i]
        # Calculationg relative velocity
        Urel  = U0*(1-a[i])/math.sin(Phi[i]) # OR PYTHAGORAS: math.sqrt((U0*(1-a[i]))**2+(Omega*dr*(i+1)*(1+aa[i]))**2)
        # Calculating drag and lift coefficiet by interpolating coefficient data
        cl = cl_data[14+math.floor(Alpha[i]*180/math.pi)]*(math.ceil(Alpha[i]*180/math.pi)-Alpha[i]*180/math.pi) + cl_data[14+math.ceil(Alpha[i]*180/math.pi)]*(Alpha[i]*180/math.pi-math.floor(Alpha[i]*180/math.pi))
        cd = cd_data[14+math.floor(Alpha[i]*180/math.pi)]*(math.ceil(Alpha[i]*180/math.pi)-Alpha[i]*180/math.pi) + cd_data[14+math.ceil(Alpha[i]*180/math.pi)]*(Alpha[i]*180/math.pi-math.floor(Alpha[i]*180/math.pi))
        # Calculation Momentum and Thrust from the forces
        dM[i] = N*0.5*Rho*Urel**2*(cl*math.sin(Phi[i])-cd*math.cos(Phi[i]))*c[i]*r[i]*dr
        dT[i] = N*0.5*Rho*Urel**2*(cl*math.cos(Phi[i])+cd*math.sin(Phi[i]))*c[i]*dr
        # Calculation of the induction factors
        Sigma = N*c[i]/(2*math.pi*r[i]) # Local solitity
        a_new[i] = 1/(1+(2*math.sin(Phi[i]))**2/(Sigma*cl*math.cos(Phi[i])))
        aa_new[i] = max(1/((4*math.cos(Phi[i])/(Sigma*cl))-1), 0)
        Difference += (a_new[i]-a[i])**2+(aa_new[i]-aa[i])**2
        a[i] = a_new[i]
        aa[i] = aa_new[i]
    evolution_a_mean.append(sum(a)/n)
    evolution_aa_mean.append(sum(aa)/n)
    print(counter, "Difference= ", Difference, "a_mean = ", evolution_a_mean[counter], "aa_mean = ", evolution_aa_mean[counter] )
    M = sum(dM)     # Calculation of total momentum
    T = sum(dT)     # Calculation of total thrust     
    P = M*Omega     # Calculation of output power
    Cp = P/Pin      # Calculation of Power coefficient
    print("Pout= ", P , "\nPin= ", Pin, "\nCp= " , Cp)
# End of the iteration
#plt.plot(evolution_a_mean)
plt.figure()
plt.plot(r, Beta/math.pi*180, label = "Beta")
plt.plot(r, Phi/math.pi*180, label = "Phi")
plt.plot(r, Alpha/math.pi*180, label = "Alpha")
plt.legend()

plt.figure()
plt.plot(r,a,label = "a")
plt.legend()

plt.figure()
plt.plot(r,aa, label = "a'")
plt.legend()

plt.figure()
plt.plot(r, dM, "g-", label='Momentum')
plt.plot(r, dT, "r-", label = 'Thrust')
plt.xlabel("Balde Element")
plt.ylabel("Force in N / Momntum in Nm")
plt.title("Momentum and thrust force on each blade eleent") 
plt.legend()
plt.show()

file = open("Alpha.txt","w") 
for num in Alpha:
    file.write(str(num) + "\n")
file.close() 
