#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e
from matplotlib.colors import BoundaryNorm

#Aufgabenteil a

#Konstanten
f = 2856e6
E = 12e6
ds = 0.01 #Schrittgröße 1cm

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

#Unterschiedliche Startbedingungen werden gesetzt
for l in range(0, 149, ):
    for n in range(0, 99):
        Ekin[0] = (l+1)*0.1* m_e*c**2
        phi[0] = - np.pi + n*np.pi/50

        for i in range(1,306):
            Ekin[i] = Ekin[i-1] + dEkin(phi[i-1])
            v[i] = V(Ekin[i])
            phi[i] = phi[i-1] + dphi(v[i])

#Endwert der kinetischen Energie wird in Abhängigkeit der Startwerte gespeichert
        Eende[l][n] = Ekin[305]


#Plotten
fig, ax = plt.subplots()

y = np.linspace(0.1, 15, 150)
x = np.linspace(-np.pi , np.pi, 100)

im = ax.pcolormesh(x, y, Eende/ e /1e6 )
fig.colorbar(im, ax=ax)
ax.set_title("kinetische Energie / MeV")

ax.set_xlabel("Phase / rad")
ax.set_ylabel("kinetische Anfangsenergie / Ruheenergie")

plt.savefig("Blatt5_Aufgabe3a_Feline,Hannah.pdf")
#weiße Felder sind Teilchen, die das Ende des Beschleunigungsmoduls nicht erreichen





