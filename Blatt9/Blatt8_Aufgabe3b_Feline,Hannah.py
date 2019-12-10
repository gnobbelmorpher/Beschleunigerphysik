#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e, mu_0
from mpl_toolkits.mplot3d import axes3d

m = 7500
b = 20**3 

x, y = np.meshgrid(np.arange(-30, 30, 1),
                      np.arange(-30, 30, 1))

u = 1e-6*m*x*y
v = 1e-6*0.5*m*(x**2-y**2)

fig = plt.figure()
ax = fig.gca()

a=np.arange(1, 31, 0.01)
print(a)

ax.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1)
ax.plot(np.sqrt(20**3/(3*a)+a**2/3), a, ) 
ax.plot(np.sqrt(20**3/(3*a)+a**2/3), -a, ) 
ax.plot(-np.sqrt(20**3/(3*a)+a**2/3), -a, ) 
ax.plot(-np.sqrt(20**3/(3*a)+a**2/3), a, ) 
ax.plot(-np.sqrt(20**3/(-3*a)+a**2/3), a, ) 
ax.plot(np.sqrt(20**3/(-3*a)+a**2/3), a, ) 
ax.plot(-np.sqrt(20**3/(-3*a)+a**2/3), -a, ) 
ax.plot(np.sqrt(20**3/(-3*a)+a**2/3), -a, ) 
ax .set_ylim([-30, 30])
ax .set_xlim([-30, 30])
#for vektor in schale:
#    ax.plot(100*vektor[0], 100*vektor[1], "o", "r")





#plt.title("Teilchenbewegung im longitudinalen Phasenraum")

#plt.xlabel(r"$\Delta \psi / rad$")
#plt.ylabel(r"$\Delta E / E$")

plt.savefig("Blatt9_Aufgabe3b_Feline,Hannah.pdf")


#stabile Trajektorien: Phase: ca. -2,5 bis 2,5
#                      Energie: bis Abweichung von 0,34 stabil
# Elektronen im instabilen Bereich verlieren immer weiter an Energie 



