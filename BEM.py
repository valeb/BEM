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

from BEM_Funcitons import LiftCoeff
from BEM_Funcitons import Performance_Wind_Turbine

# Calling the functions ###############################################################

Airfoilnames = ['./2DAirfoilDataFiles/NACA 0615A.txt','./2DAirfoilDataFiles/NACA 1412.txt','./2DAirfoilDataFiles/NACA 23021.txt','./2DAirfoilDataFiles/NACA 63012A.txt','./2DAirfoilDataFiles/NACA 64008A.txt', './2DAirfoilDataFiles/DU 06W200.txt']

n = 100

Mmax = [0]*n
Foilmax = [len(Airfoilnames)]*n
for i  in range(len(Airfoilnames)) :
    M = Performance_Wind_Turbine(10, Airfoilnames[i], False)
    print(M)
    for j in range(n):
        if M[j] > Mmax[j] :
            Mmax[j] = M[j]
            Foilmax[j] = i

print(Foilmax)

for foil in Foilmax :
    if foil < len(Airfoilnames) :    
        print(Airfoilnames[foil])      
    else :
        print('Cylinder')
plt.show()
