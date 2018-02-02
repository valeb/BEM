#########################################################################
#                                                                       #
#  Performance analysis of a wind turbine for a given incoming wind     #
#                                                                       #
#  Alija Bajramovic                                                     #
#  Valentin Bernard                                                     #
#  Naven Goutham                                                        #
#                                                                       #
#########################################################################

# Lbraries #############################################################################

import math
import matplotlib.pyplot as plt
import numpy as np
import tkinter.filedialog

# Functions ############################################################################

def LiftCoeff (alpha):
    if alpha<21:
        Cl = 0.42625 + 0.11628 * alpha - 0.00063973 * alpha**2 - 8.712 * 10**(-5) * alpha**3 - 4.2576 * 10**(-6)*alpha**4
    else:
        Cl = 0.95
    return Cl

# Definition of Variales ###############################################################

# Physical Parameters
U0 = 10                             # Incoming wind speed in m/s
Rho = 1.225                         # Air density in kg/m^3
R = 50                              # length of the blade in m
N = 3                               # Number of blades
Lambda = 8                          # Tip speed ratio

# Computational Parmeters
n = 100                             # Number of blade elements
dr = R/100                          # Blade element length
r_e = np.arange(dr, R+dr, dr)       # End point radius of every element
r = r_e -(dr/2)                     # Midpoint radius
Lambda_r = r*Lambda/R               # Local tip speed ratio

# Blade characteristics
c = [1] * n                         # Airfoil chord in m as function of r
Beta = np.genfromtxt("Beta.txt")    # Twist angle in radiants as function of r

# Derived quantities
Omega = Lambda*U0/R                 # Rotor angular velocity in radiants/s
Area = R**2*math.pi                 # Rotor area in m^2
Pin = 1/2*U0**3*Rho*Area            # Incoming wind power in W
Sigma = np.zeros(n)                 # Local solitity
for i in range(n):
    Sigma[i] = min(N*c[i]/(2*math.pi*r[i]), 1)


# Results
Cp = 0.0                            # Power Coefficient
Ct = 0.0                            # Thrust coefficient
dL = [0.0] * n                      # Lift on each elment in N
dD = [0.0] * n                      # Drag on each elment in N
dM = [0.0] * n                      # Momentm on each elment in Nm
M = 0.0                             # Total momentum in Nm
P = 0.0                             # Power in W
dT = [0.0] * n                      # Thrust force on each element in N
T = 0.0                             # Total thrust force
Phi = np.zeros(n)                   # Angle of incoming wind in radiants
Alpha = np.zeros(n)                 # Angle of attack in radiants

# Main Program #######################################################################

a = [1/3]*n                         # Linear induction factor, initial guess
aa = [0.0]*n                        # Angular induction factor, initial guess

# Variables for the iteration
Difference = 999
counter = 0

# Loop for the iteration
while (Difference > 10**(-10)) :
    Difference = 0.0
    counter += 1
    for i in range(n):
        # Calculationg the angle of the incoming wind
        Phi[i] = math.atan((1-a[i])/(1+aa[i])/Lambda_r[i])
        Alpha[i] = Phi[i]-Beta[i]
        # Calculationg relative velocity
        Urel  = U0*(1-a[i])/math.sin(Phi[i])
        # Calculating drag and lift coefficiet from empiric equation
        Cl = LiftCoeff(Alpha[i])
        cd = 0
        # Calculation Momentum and Thrust from the forces
        dM[i] = N*0.5*Rho*Urel**2*(Cl*math.sin(Phi[i])-cd*math.cos(Phi[i]))*c[i]*r[i]*dr
        dT[i] = N*0.5*Rho*Urel**2*(Cl*math.cos(Phi[i])+cd*math.sin(Phi[i]))*c[i]*dr
        # Calculation of the induction factors
        a_new = 1/(1+(2*math.sin(Phi[i]))**2/(Sigma[i]*Cl*math.cos(Phi[i])))
        aa_new = max(0, 1/((4*math.cos(Phi[i])/(Sigma[i]*Cl))-1))
        if i > 5: Difference += (a_new-a[i])**2+(aa_new-aa[i])**2
        a[i] = a_new
        aa[i] = aa_new
    M = sum(dM)     # Calculation of total momentum
    T = sum(dT)     # Calculation of total thrust     
    P = M*Omega     # Calculation of output power
    Cp = P/Pin      # Calculation of Power coefficient
    print(counter, ":", "Pout= ", P , "Pin= ", Pin, "Cp= " , Cp)
# End of the iteration

# Plots ################################################################################

# Plot Angles
plt.figure()
plt.plot(r, Beta/math.pi*180, label = "Beta")
plt.plot(r, Phi/math.pi*180, label = "Phi")
plt.plot(r, Alpha/math.pi*180, label = "Alpha")
plt.legend()
#plt.xlim(left = 3)

# Plot the 2 induction factors along the blade

# Two subplots, unpack the axes array immediately
f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
ax1.plot(r,a)
ax2.plot(r,aa)
ax1.set_title("a(r)")
ax2.set_title("a'(r)")

# Plot Torque and Thrust on the blade elements
plt.figure()
plt.plot(r, dM, "g-", label='Momentum')
plt.plot(r, dT, "r-", label = 'Thrust')
plt.xlabel("Balde Element")
plt.ylabel("Force in N / Momntum in Nm")
plt.title("Momentum and thrust force on each blade eleent") 
plt.legend()
#plt.xlim(left = 3)

# File output ############################################################################

# Export resulting Alpha along the blade
file = open("Alpha.txt","w") 
for num in Alpha:
    file.write(str(num) + "\n")
file.close()


plt.show() 
