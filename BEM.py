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

Cp = [0]*300
plt.figure
for i in range(40,300):
    Cp[i] = Performance_Wind_Turbine(i/10, False)
    plt.plot(i/10, Cp[i], 'b*')

#print(Performance_Wind_Turbine(10, True))

plt.show()
