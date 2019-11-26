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
W = 570e3 #eV
phis = np.pi - np.arcsin(W/(V))



def deltaphi(dE, a0, a1, a2):
    return omega * (a0 + a1 * dE + a2 * (dE)**2)

def deltaE(dphi):
    return e * V / (T * e*E) * ( np.sin(phis + dphi) - np.sin(phis) )


dE = np.zeros(500)
dphi = np.zeros(500)

a0 = 0
a1 = 8e-5
a2 = 1e-3

for n in range(1, 20):
    for k in range(0,9):
    #k= (n-1)*1001
        dE[0] = 0.002*n
        dphi[0] = k * 2 * np.pi /10 - np.pi
        for i in range(1, 500):

            j1 = T * deltaphi(dE[i-1], a0, a1, a2)
            k1 = T * deltaE(dphi[i-1])
            j2 = T * deltaphi(dE[i-1]+k1/2, a0, a1, a2)
            k2 = T * deltaE(dphi[i-1]+j1/2)
            j3 = T * deltaphi(dE[i-1]+k2/2, a0, a1, a2)
            k3 = T * deltaE(dphi[i-1]+j2/2)
            j4 = T * deltaphi(dE[i-1]+k3, a0, a1, a2)
            k4 = T * deltaE(dphi[i-1]+j3)

            dphi[i] = dphi[i-1] + 1/6 * (j1+ 2*j2 + 2*j3 + j4)
            dphi[i] = dphi[i]%(2*np.pi) #Wertebereich von 0 bis 2pi
            if dphi[i] > np.pi:
                dphi[i] = dphi[i] - 2*np.pi #wertebereich von -pi bis pi
            dE[i] = dE[i-1] + 1/6 * (k1+ 2*k2 + 2*k3 + k4)
        plt.plot(dphi, dE, ",")



plt.title(r"$\alpha_2 = 10^{-3}, \alpha_1 = 8 \cdot 10^{-5},  $")

plt.xlabel(r"$\Delta \psi / rad$")
plt.ylabel(r"$\Delta E / E$")

plt.savefig("Blatt7_Aufgabe3b10_Feline,Hannah.pdf")


#stabile Trajektorien: Phase: ca. -2,5 bis 2,5
#                      Energie: bis Abweichung von 0,34 stabil
# Elektronen im instabilen Bereich verlieren immer weiter an Energie 



