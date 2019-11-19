#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e
from matplotlib.colors import BoundaryNorm

#Aufgabenteil a

#Konstanten
E = 1.7e9 #eV
V = 1.2e6 #V
T = 800e-9 #s
f = 500e6 #Hz
omega = 2* np.pi * f
a = 8e-4 #alpha
W = 170e3 #eV
phis = np.pi - np.arcsin(W/(V))
print(phis)


#Änderungen

def deltaphi(dE):
    return omega * a * dE * T

def deltaE(dphi):
    return e * V / (T * e*E) * ( np.sin(phis + dphi) - np.sin(phis) ) * T


#Arrays erstellen
dE = np.zeros(1000)
dphi = np.zeros(1000)
dphi[0] = 0


#Unterschiedliche Startbedingungen werden gesetzt
for n in range(1, 20):
    #k= (n-1)*1001
    dE[0] = 0.002*n
    for i in range(1, 1000):
        dphi[i] = dphi[i-1] + deltaphi(dE[i-1])
        dphi[i] = dphi[i]%(2*np.pi)
        if dphi[i] > np.pi:
            dphi[i] = dphi[i] - 2*np.pi
        dE[i] = dE[i-1] + deltaE(dphi[i])
    plt.plot(dphi, dE, ".")


#Plotten
#fig, ax = plt.subplots()
#plt.plot(dphi, dE)
#y = np.linspace(0.1, 15, 150)
#x = np.linspace(-np.pi , np.pi, 100)

#im = ax.pcolormesh(x, y, Eende/ e /1e6 )
#fig.colorbar(im, ax=ax)
#ax.set_title("kinetische Energie / MeV")

#ax.set_xlabel("Phase / rad")
#ax.set_ylabel("kinetische Anfangsenergie / Ruheenergie")

plt.savefig("Blatt5_Aufgabe3a_Feline,Hannah.pdf")
#weiße Felder sind Teilchen, die das Ende des Beschleunigungsmoduls nicht erreichen





