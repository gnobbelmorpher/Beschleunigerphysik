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

mass =  np.sqrt( 2 * (511e-6)**2 + 2 * np.outer(energy1, energy2) + 2 * np.outer(np.sqrt(energy1**2 - (511e-6)**2), np.sqrt(energy2**2 - (511e-6)**2)))

#Aufgabenteil b

plt.plot(energy1, mass[0])
plt.xlabel(r"Teilchenenergie / GeV")
plt.ylabel(r"invariante Masse / GeV/c²")
plt.savefig("Blatt2_Aufgabe3b_Feline,HannahKorrektur.pdf")

#man sieht den parabolischen verlauf der invarianten masse

#Aufgabenteil c


fig = plt.figure()
ax = Axes3D(fig)
x, y = np.meshgrid(energy1, energy2,)
ax.plot_surface(x, y, mass, rstride=1, cstride=1, cmap=cm.viridis)
plt.xlabel(r"Energie Teilchen 1 / GeV")
plt.ylabel(r"Energie Teilchen 2 / GeV")
ax.set_zlabel(r"invariante Masse / GeV/c²")
plt.savefig("Blatt2_Aufgabe3c_Feline,HannahKorrektur.pdf")

fig, ax = plt.subplots()
im = ax.imshow(mass)

cbar = ax.figure.colorbar(im, ax=ax, )
cbar.ax.set_ylabel("Invariante Masse / GeV/c²", rotation=-90, va="bottom")

a = np.linspace(0, 99, 7 )
b = np.linspace(0, 6, 7)

# We want to show all ticks...
ax.set_xticks(a)
ax.set_yticks(a)
# ... and label them with the respective list entries
ax.set_xticklabels(b)
ax.set_yticklabels(b)
plt.xlabel(r"Energie Teilchen 1 / GeV")
plt.ylabel(r"Energie Teilchen 2 / GeV")
plt.savefig("Blatt2_Aufgabe3c_Feline,HannahHeatmap.pdf")


#auch hier ein parabolischer verlauf der invarianten masse






