#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e, mu_0
from mpl_toolkits.mplot3d import axes3d


def biotsavart(I, x):
    return - mu_0 /(2 * np.pi) * np.cross(I,x / np.dot(x,x))

In = [0,0,2]
Is = [0,0,-2]
schale=np.zeros((40 ,3))
r = 0.04
phi0 = np.pi * 0.3

#platzieren der drähte in Polarkoordinaten
for i in range(0,20):
    schale[i] = [r, -phi0+i*phi0/10+np.pi/2, 0]

for i in range(20,40):
    schale[i] = [r, np.pi-phi0+(i-20)*phi0/10+np.pi/2, 0]

#umrechnen in kartesische koordinaten
for vektor in schale:
    r = vektor[0] 
    phi = vektor[1]
    vektor[0] = r * np.sin(phi)
    vektor[1] = r * np.cos(phi)

B = np.zeros((100,100,3))

for i in range(0,20):
    for k in range(1,99):
        for n in range(1,99): 
            x = [0.001*k -0.05 - schale[i,0], 0.001*n - 0.05- schale[i,1], 0]
            B[k,n] = B[k,n] + biotsavart(Is, x)

for i in range(20,40):
    for k in range(0,100):
        for n in range(0,100): 
            x = [0.001*k - 0.05 - schale[i,0], 0.001*n - 0.05 - schale[i,1], 0]
            B[k,n] = B[k,n] + biotsavart(In, x)


fig = plt.figure()
ax = fig.gca()

# Make the grid
x, y = np.meshgrid(np.arange(-50, 50, 1),
                      np.arange(-50, 50, 1))

u = 1e4*B[x-50, y-50, 0]
v = 1e4*B[x-50, y-50, 1]

ax.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1)
ax.plot(1000*schale[:,0], 1000*schale[:,1], 'x') 
#for vektor in schale:
#    ax.plot(100*vektor[0], 100*vektor[1], "o", "r")





#plt.title("Teilchenbewegung im longitudinalen Phasenraum")

#plt.xlabel(r"$\Delta \psi / rad$")
#plt.ylabel(r"$\Delta E / E$")

plt.savefig("Blatt8_Aufgabe3a_Feline,Hannah.pdf")


#stabile Trajektorien: Phase: ca. -2,5 bis 2,5
#                      Energie: bis Abweichung von 0,34 stabil
# Elektronen im instabilen Bereich verlieren immer weiter an Energie 



