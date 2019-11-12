#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e

#Aufgabenteil a

#Konstanten
f = 2856e6
E = 12e6
ds = 0.01

#Änderungen

def dphi(v):
    return 2 * np.pi * f / c * (v - c) * ds / v

def dEkin(phi):
    return e * E * np.cos(phi) * ds

def V(Ekin):
    return c* np.sqrt(1-(m_e*c**2)/(Ekin + m_e*c**2))

#Arrays erstellen
Ekin = np.zeros(306)
v = np.zeros(306)
phi = np.zeros(306)
Eende = np.zeros((150,100))

for l in range(0, 149, ):
    for n in range(0, 99):
        Ekin[0] = (l+1)*0.1* m_e*c**2
        phi[0] = - np.pi + n*np.pi/50

        for i in range(1,306):
            Ekin[i] = Ekin[i-1] + dEkin(phi[i-1])
            v[i] = V(Ekin[i])
            phi[i] = phi[i-1] + dphi(v[i])

        Eende[l][n] = Ekin[305]



#print(Eende)
fig, ax = plt.subplots()
im = ax.imshow(Eende / e / 1e6)

cbar = ax.figure.colorbar(im, ax=ax, )
cbar.ax.set_ylabel("kinetische Energie / MeV", rotation=-90, va="bottom")

plt.savefig("Blatt5_Aufgabe3a_Feline,Hannah.pdf")





