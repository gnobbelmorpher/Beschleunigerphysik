#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e, mu_0
from mpl_toolkits.mplot3d import axes3d


def biotsavart(I, x):
    return - mu_0 /(2 * np.pi) * np.cross(I,x / np.dot(x,x))

In = [0,0,1]
Is = [0,0,-1]
schale=np.zeros((40 ,3))
r = 0.04
phi0 = np.pi * 0.5

#platzieren der drähte in Polarkoordinaten
for i in range(0,19):
    schale[i] = [r, -phi0+i*phi0/10, 0]

for i in range(20,39):
    schale[i] = [r, np.pi-phi0+i*phi0/10, 0]

#umrechnen in kartesische koordinaten
for vektor in schale:
    r = vektor[0] 
    phi = vektor[1]
    vektor[0] = r * np.sin(phi)
    vektor[1] = r * np.cos(phi)
    print(vektor)

B = np.zeros((100,100,3))

#südpol
for i in range(0,19):
    for k in range(1,99):
        for n in range(1,99): 
            x = [0.01*k - schale[i,0], 0.01*n - schale[i,1], 0]
            B[k,n] = B[k,n] + biotsavart(Is, x)

#nordpol
for i in range(20,39):
    for k in range(1,99):
        for n in range(1,99): 
            x = [0.001*k - schale[i,0], 0.001*n - schale[i,1], 0]
            B[k,n] = B[k,n] + biotsavart(In, x)


fig = plt.figure()
ax = fig.gca()

# Make the grid
x, y = np.meshgrid(np.arange(-50, 50, 1),
                      np.arange(-50, 50, 1))

u = 1e4*B[x, y, 0]
v = 1e4*B[x, y, 1]

ax.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1)
for vektor in schale:
    plt.plot(100*vektor[0], 100*vektor[1])





#plt.title("Teilchenbewegung im longitudinalen Phasenraum")

#plt.xlabel(r"$\Delta \psi / rad$")
#plt.ylabel(r"$\Delta E / E$")

plt.savefig("Blatt8_Aufgabe3_Feline,Hannah.pdf")


#stabile Trajektorien: Phase: ca. -2,5 bis 2,5
#                      Energie: bis Abweichung von 0,34 stabil
# Elektronen im instabilen Bereich verlieren immer weiter an Energie 



