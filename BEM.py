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

Cp = [0]*30
plt.figure
#for i in range(30):
Cp[i] = Performance_Wind_Turbine(10, True)#i+1, False)
plt.plot(i, Cp[i], 'b*')

plt.show()
