#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

#Aufgabeinteil a

energy1 = np.linspace(511e-6, 6, 100 )
energy2 = np.linspace(511e-6, 6, 100 )

#Berechnung der invarianten Masse

mass =  np.sqrt( 2 * (511e-6)**2 + 2 * np.outer(energy1, energy2) + np.outer(np.sqrt(energy1**2 - (511e-6)**2), np.sqrt(energy2**2 - (511e-6)**2)))
#mass =  np.sqrt(np.outer(energy1, energy2)- (np.sqrt(energy1**2 )))
#mass =  np.sqrt(np.outer(energy1, energy2))

#Aufgabenteil b

plt.plot(energy1, mass[0])
plt.xlabel(r"Teilchenenergie / GeV")
plt.ylabel(r"invariante Masse / GeV/c²")
plt.savefig("Blatt2_Aufgabe3b_Feline,Hannah.pdf")

#Aufgabenteil c


fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(energy1, energy2, mass, rstride=1, cstride=1, cmap=cm.viridis)
plt.xlabel(r"Energie Teilchen 1 / GeV")
plt.ylabel(r"Energie Teilchen 2 / GeV")
ax.set_zlabel(r"invariante Masse / GeV/c²")
plt.savefig("Blatt2_Aufgabe3c_Feline,Hannah.pdf")






