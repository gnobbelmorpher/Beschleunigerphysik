#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e, mu_0
from mpl_toolkits.mplot3d import axes3d

g = 75

x, y = np.meshgrid(np.arange(-30, 30, 1),
                      np.arange(-30, 30, 1))

u = 1e-3*g*y
v = 1e-3*g*x

fig = plt.figure()
ax = fig.gca()

a=np.arange(1, 31, 0.1)
print(a)

ax.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1)
ax.plot(a, 200/a, ) 
ax.plot(-a, 200/a, ) 
ax.plot(a, -200/a, ) 
ax.plot(-a, -200/a, ) 
ax .set_ylim([-30, 30])
#for vektor in schale:
#    ax.plot(100*vektor[0], 100*vektor[1], "o", "r")





#plt.title("Teilchenbewegung im longitudinalen Phasenraum")

#plt.xlabel(r"$\Delta \psi / rad$")
#plt.ylabel(r"$\Delta E / E$")

plt.savefig("Blatt9_Aufgabe3a_Feline,Hannah.pdf")


#stabile Trajektorien: Phase: ca. -2,5 bis 2,5
#                      Energie: bis Abweichung von 0,34 stabil
# Elektronen im instabilen Bereich verlieren immer weiter an Energie 



