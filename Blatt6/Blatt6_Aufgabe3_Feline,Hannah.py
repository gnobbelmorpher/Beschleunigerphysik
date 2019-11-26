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

#leider ist die Bedennung etwas doof: deltaphi bei mir ist d(\Dela \psi)
#und dphi ist \delta \phi (bei E identisch)


def deltaphi(dE):
    return omega * a * dE * T

def deltaE(dphi):
    return e * V / (T * e*E) * ( np.sin(phis + dphi) - np.sin(phis) ) * T


dE = np.zeros(250)
dphi = np.zeros(250)
dphi[0] = 0


for n in range(1, 20):
    #k= (n-1)*1001
    dE[0] = 0.002*n
    for i in range(1, 250):
        dphi[i] = dphi[i-1] + deltaphi(dE[i-1])
        dphi[i] = dphi[i]%(2*np.pi) #Wertebereich von 0 bis 2pi
        if dphi[i] > np.pi:
            dphi[i] = dphi[i] - 2*np.pi #wertebereich von -pi bis pi
        dE[i] = dE[i-1] + deltaE(dphi[i])
    plt.plot(dphi, dE, ".")



plt.title("Teilchenbewegung im longitudinalen Phasenraum")

plt.xlabel(r"$\Delta \psi / rad$")
plt.ylabel(r"$\Delta E / E$")

plt.savefig("Blatt5_Aufgabe3a_Feline,Hannah.pdf")


#stabile Trajektorien: Phase: ca. -2,5 bis 2,5
#                      Energie: bis Abweichung von 0,34 stabil
# Elektronen im instabilen Bereich verlieren immer weiter an Energie 



