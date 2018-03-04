
#########################################################################
#                                                                       #
#  Oprimization of a wind turbine for a given incoming wind             #
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
from os import listdir
from os.path import isfile, join

from BEM_Funcitons import Performance_Wind_Turbine

# Global Variables ####################################################################

# Physical Parameters
U0 = 10                             # Incoming wind speed in m/s
Rho = 1.225                         # Air density in kg/m^3
R = 50                              # length of the blade in m
N = 3                               # Number of blades
Omega = 1.6                         # Rotor angular velocity to obtain Lambda = 8 at 10m/s
Area = R**2*math.pi                 # Rotor area in m^2
Lambda = R*Omega/U0
Pin = 1/2*U0**3*Rho*Area            # Incoming wind power in W

# Computational Parmeters
n = 100                             # Number of blade elements
dr = R/n                      # Blade element length
r = np.arange(dr, R+dr, dr)     # Mid point radius of every element

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

def Iterate_a_aa (AlphaData, ClData, Lambda_r, Beta, Sigma) :
    # Variables for the iteration
    Difference = 10
    a = 1/3
    aa = 0
    counter = 0;
    Limit = 10**(-4)
    # Loop for the iteration
    while Difference > Limit :
    # Calculationg the angle of the incoming wind
        Phi = math.atan((1-a)/(1+aa)/Lambda_r)
        Alpha = Phi*180/math.pi - Beta
    # Calculating drag and lift coefficiet from empiric equation
        Cl = Coefficient(Alpha, AlphaData, ClData)
        Cd = Coefficient(Alpha, AlphaData, CdData)
        if Cl == 0 : Cl += 0.0001
    # Calculation of the new induction factors
        a_new = 1/(1+(2*math.sin(Phi))**2/(Sigma*Cl*math.cos(Phi)))
        aa_new = 1/((4*math.cos(Phi)/(Sigma*Cl))-1)
        Difference = (a_new-a)**2+(aa_new-aa)**2
        a = a_new
        aa = aa_new
        counter += 1
        if counter > 50 : return([9,9,9,9,9]) #raise ValueError('a aa did not converge, Beta =', Beta, 'Sigma =', Sigma) 
    # End of the iteration 
    return([a, aa, Phi, Cl, Cd] )

# Calling the functions ###############################################################

Airfoils = [join('./2DAirfoilDataFiles',f) for f in listdir('./2DAirfoilDataFiles') if isfile(join('./2DAirfoilDataFiles', f)) and f[-1] == 't']
#Airfoils = [join('./2DAirfoilDataFiles/NACA',f) for f in listdir('./2DAirfoilDataFiles/NACA') if isfile(join('./2DAirfoilDataFiles/NACA', f)) and f[-1] == 't']

print(Airfoils)


# Define best performance variables!
Mmax = [0]*n
Cmax = [2]*n
BetaMax = [0]*n
Fmax = [99]*n

for j in range(1,len(Airfoils)) :
    print(j, Airfoils[j])
    AlphaData = np.genfromtxt(Airfoils[j], usecols=0)
    ClData = np.genfromtxt(Airfoils[j], usecols=1)
    CdData = np.genfromtxt(Airfoils[j], usecols=2)
    for i in range(n-1, 4, -1) :
        r = (i+1)*dr                # Local radius
        Lambda_r = r*Lambda/R       # Local tip speed ratio
        for Beta in np.arange(-20, 20, 1) :
            for c in np.arange(0.1, min(math.ceil((2*math.pi*r)/N), 5), 1) :
                Sigma = N*c/(2*math.pi*r)
                [a, aa, Phi, Cl, Cd] = Iterate_a_aa(AlphaData, ClData, Lambda_r, Beta, Sigma)
                if a == 9 : break
                # Calculationg relative velocity
                Urel  = U0*(1-a)/math.sin(Phi)
                # Calculation Momentum from the forces
                M = N*0.5*Rho*Urel**2*(Cl*math.sin(Phi)-Cd*math.cos(Phi))*c*r*dr
                if M > Mmax[i] :
                    Mmax[i] = M
                    Cmax[i] = c
                    BetaMax[i] = Beta
                    Fmax[i] = j
        for Beta in np.arange(BetaMax[i]-2, BetaMax[i]+2, 0.1) :
            for c in np.arange(Cmax[i]-2, Cmax[i]+2, 0.1) :
                Sigma = N*c/(2*math.pi*r)
                [a, aa, Phi, Cl, Cd] = Iterate_a_aa(AlphaData, ClData, Lambda_r, Beta, Sigma)
                if a == 9 : break
                # Calculationg relative velocity
                Urel  = U0*(1-a)/math.sin(Phi)
                # Calculation Momentum and Thrust from the forces
                M = N*0.5*Rho*Urel**2*(Cl*math.sin(Phi)-Cd*math.cos(Phi))*c*r*dr
                if M >= Mmax[i] :
                    Mmax[i] = M
                    Cmax[i] = c
                    BetaMax[i] = Beta
                    Fmax[i] = j
        for Beta in np.arange(BetaMax[i]-0.2, BetaMax[i]+0.2, 0.01) :
            for c in np.arange(Cmax[i]-0.2, Cmax[i]+0.2, 0.01) :
                Sigma = N*c/(2*math.pi*r)
                [a, aa, Phi, Cl, Cd] = Iterate_a_aa(AlphaData, ClData, Lambda_r, Beta, Sigma)
                if a == 9 : break
                # Calculationg relative velocity
                Urel  = U0*(1-a)/math.sin(Phi)
                # Calculation Momentum and Thrust from the forces
                M = N*0.5*Rho*Urel**2*(Cl*math.sin(Phi)-Cd*math.cos(Phi))*c*r*dr
                if M >= Mmax[i] :
                    Mmax[i] = M
                    Cmax[i] = c
                    BetaMax[i] = Beta
                    Fmax[i] = j


P = sum(Mmax)*Omega
Cp = P/Pin

# Outputs

for i in range(5,n) :
    print(i, 'c =', Cmax[i], 'Beta =', BetaMax[i], 'Foil =', Airfoils[Fmax[i]][21:-4],'(',Fmax[i],')', )

print('cp =', Cp)


plt.figure()
plt.plot(Cmax, label = "c")
plt.legend()

plt.figure()
plt.plot(BetaMax, label = "Beta")
plt.legend()

# File output ####################################################################

# Export resulting c and Beta
file = open("Beta.txt","w") 
for num in BetaMax:
    file.write(str(num) + "\n")
file.close()

file = open("c.txt","w") 
for num in Cmax:
    file.write(str(num) + "\n")
file.close()

file = open("Foils.txt","w") 
for num in Fmax:
    if num == 99 :
        file.write( 'Cylinder.txt' + "\n")
    else :
        file.write(Airfoils[num] + "\n")
file.close()

plt.show()
