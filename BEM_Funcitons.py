###########################################################################
#                                                                         #
#  Performance analysis of a wind turbine for a given incoming wind       #
#  and angular velocity using blade element momentum theory               #
#                                                                         #
#  Alija Bajramovic                                                       #
#  Valentin Bernard                                                       #
#  Naveen Goutham                                                         #
#                                                                         #
###########################################################################

import math
import matplotlib.pyplot as plt
import numpy as np
import tkinter.filedialog


# Functions ############################################################################

def Coefficient (alpha, AlphaData, CoeffData): 
    alpha = alpha*180/math.pi
    # Interpolate data to obtain Cl or Cd for given alpha in degrees
    i = 0
    if alpha < AlphaData[0]: return CoeffData[0]
    if alpha > AlphaData[-1]: return CoeffData[-1]
    while not (AlphaData[i] <= alpha and AlphaData[i+1] > alpha):
        i += 1
    coeff = (CoeffData[i+1]*(alpha-AlphaData[i]) + CoeffData[i]*(AlphaData[i+1]-alpha))/(AlphaData[i+1]-AlphaData[i]) 
    return coeff

def LiftCoeff (alpha, AlphaData, ClData):
    alpha = alpha*180/math.pi
    if alpha<21:
        Cl = 0.42625 + 0.11628 * alpha - 0.00063973 * alpha**2 - 8.712 * 10**(-5) * alpha**3 - 4.2576 * 10**(-6)*alpha**4
    else:
        Cl = 0.95
    return Cl

def Performance_Wind_Turbine(Incoming_Wind, foilfile, printing):
    # Evaluate performance of a wind turbine given the incoming wind.
    # The wind turbine is defined by files Beta.txt, c.txt and foilfile!

    # Definition of Variables       ##########################################################

    # Physical Parameters
    #U0 = 10                             # Incoming wind speed in m/s
    U0 = Incoming_Wind
    Rho = 1.225                         # Air density in kg/m^3
    R = 50                              # length of the blade in m
    R_i = 5                             # Inner radius in m
    N = 3                               # Number of blades
    Omega = 1.6                         # Rotor angular velocity optimized for Lambda = 8

    # Computational Parmeters
    n = 100                             # Number of blade elements
    dr = (R-R_i)/n                      # Blade element length
    r = np.arange(dr+R_i, R+dr, dr)     # Mid point radius of every element
    Lambda = R*Omega/U0
    Lambda_r = r*Lambda/R               # Local tip speed ratio
    
    # Blade characteristics
    c = np.genfromtxt("c.txt")          # Airfoil chord in m as function of r
    Beta = np.genfromtxt("Beta.txt")    # Twist angle in radians as function of r
    AlphaData = np.genfromtxt(foilfile, usecols=0)
    ClData = np.genfromtxt(foilfile, usecols=1)
    CdData = np.genfromtxt(foilfile, usecols=2)
    
    # Derived quantities
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
    for i in range(10,n):
        # Variables for the iteration
        Difference = 999
        # Loop for the iteration
        while (Difference > 10**(-5)) :
        # Calculationg the angle of the incoming wind
            Phi[i] = math.atan((1-a[i])/(1+aa[i])/Lambda_r[i])
            Alpha[i] = Phi[i]-Beta[i]
        # Calculating drag and lift coefficiet from empiric equation
            Cl = Coefficient(Alpha[i], AlphaData, ClData)
            #Cd = Coefficient(Alpha[i], AlphaData, CdData)
        # Calculation of the new induction factors
            a_new = 1/(1+(2*math.sin(Phi[i]))**2/(Sigma[i]*Cl*math.cos(Phi[i])))
            aa_new = 1/((4*math.cos(Phi[i])/(Sigma[i]*Cl))-1)
            Difference = (a_new-a[i])**2+(aa_new-aa[i])**2
            a[i] = a_new
            aa[i] = aa_new
        # End of the iteration, a and aa are now final for one element

        # End of the iteration, a and aa are now final for one element
    # Iteration finished for every blade element ###########################################

    # Calculate the Forces and power of the turbine
    for i in range(10,n):
    # Calculationg relative velocity
        Urel  = U0*(1-a[i])/math.sin(Phi[i])
    # Calculating drag and lift coefficiet from empiric equation
        Cl = Coefficient(Alpha[i], AlphaData, ClData)
        Cd = Coefficient(Alpha[i], AlphaData, CdData)
    # Calculation Momentum and Thrust from the forces
        dM[i] = N*0.5*Rho*Urel**2*(Cl*math.sin(Phi[i])-Cd*math.cos(Phi[i]))*c[i]*r[i]*dr
        dT[i] = N*0.5*Rho*Urel**2*(Cl*math.cos(Phi[i])+Cd*math.sin(Phi[i]))*c[i]*dr
    M = sum(dM)        # Calculation of total momentum
    T = sum(dT)        # Calculation of total thrust
    P = M*Omega        # Calculation of output power
    Cp = P/Pin         # Calculation of Power coefficient


    print("Pout= ", P , "Pin= ", Pin, "Cp= " , Cp)

    # Plots ############################################################################

    if printing:
        
        # Plot Angles
        plt.figure()
        plt.plot(r, Beta/math.pi*180, label = "Beta")
        plt.plot(r, Phi/math.pi*180, label = "Phi")
        plt.plot(r, Alpha/math.pi*180, label = "Alpha")
        plt.legend()
        #plt.xlim(left = 3)

        # Plot Cl Curve
        plt.figure()
        plt.plot(AlphaData, ClData, label = "Cl(alpha)")
        plt.plot(AlphaData, CdData, label = "Cd(alpha)")
        plt.legend()


        # Plot the 2 induction factors along the blade

        # Two subplots, unpack the axes array immediately
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
        ax1.plot(r,a)
        ax2.plot(r,aa)
        ax1.set_title("a(r)")
        ax2.set_title("a'(r)")

        # File output ####################################################################

        # Export resulting Alpha along the blade
        file = open("a.txt","w") 
        for num in a:
            file.write(str(num) + "\n")
        file.close()

        file = open("aa.txt","w") 
        for num in aa:
            file.write(str(num) + "\n")
        file.close()

        file = open("Alpha.txt","w") 
        for num in Alpha:
            file.write(str(num) + "\n")
        file.close()

        plt.show()
    return dM


