#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

#Aufgabenteil a
#nicht-relativistisch
f = 20 * 10**6
q = 1.6 * 10**(-19) 
m = 1.7 * 10**(-27)
n = np.linspace(0,400,401)

#B-Feld berechnen
B = f * 2 * np.pi *  m / q

#Dauer halber Umlauf
T = np.pi * m / (q * B) 

#Zeitpunkt des Durchflugs durch das E-Feld
t = T  * n

# kinetische Enerie in keV
Ekin =  n * 50 * np.absolute(np.cos(2 * np.pi* f * t)) + 50


#aktueller Bahnradius
r = np.sqrt(Ekin * 1e3 * q / (2 * m)) / (np.pi * f)

B
print(n)
print(np.absolute(np.cos( np.pi* n)))
print(Ekin)
plt.subplot(311)
plt.xlabel("Anzahl halbe Umläufe")
plt.ylabel("Kinetische Energie / keV")
plt.plot(n, Ekin)

plt.subplot(312)
plt.xlabel("Anzahl halbe Umläufe")
plt.ylabel("Dauer halber Umlauf / s")
plt.plot(n, np.linspace(T, T, 401))

plt.subplot(313)
plt.xlabel("Anzahl halbe Umläufe")
plt.ylabel("Radius / m")
plt.plot(n, r)

plt.savefig("Blatt4_Aufgabe3a_Feline,Hannah.pdf")
plt.clf()

#Aufgabenteil b
#relativistisch

#Ruheenergie- und Masse
#Energie in keV
E0 = 938.27 * 10**3
m0 = 1.7 * 10**(-27)

#array der richtigen Größe wird erzeugt
E = np.copy(n)
m = np.copy(n)
T = np.copy(n)
g = np.copy(n)
Ekin = np.copy(n)
r = np.copy(n)
t = np.copy(n)

#Anfangwerte werden gesetzt
E[0] = E[-1] = 50 + E0

t[0] = 0

for i in range(401):
#Bestimmung Gammafaktor und "aktuelle Masse"
#hier gehts kaputt g ist viel zu klein
    g[i] = E[i-1]/E0
    m[i] = m0 * g[i]

#Zeit halber Umlauf und Zeitpunkt des Durchlaufs
    T[i] = np.pi * m[i] / (q * B)
    t[i] = t[i-1] + T[i]
    t[0]= 0

#Berechnung aktueller Gesamtenergie
    E[i] = E[i-1] + 50 * np.absolute(np.cos(2 * np.pi * f* t[i]))
    E[0] = 50 + E0

# kinetische Enerie in keV 
    Ekin[i] = E[i]-E0

#aktueller Bahnradius

    r[i] = np.sqrt( Ekin[i] *10**-3 * q/(2*m[i])) *2* T[i] /  np.pi 

plt.subplot(311)
plt.xlabel("Anzahl halbe Umläufe")
plt.ylabel("Kinetische Energie / keV")
plt.plot(n[1:-1], Ekin[1:-1])

plt.subplot(312)
plt.xlabel("Anzahl halbe Umläufe")
plt.ylabel("Dauer halber Umlauf / s")
plt.plot(n, T)

plt.subplot(313)
plt.xlabel("Anzahl halbe Umläufe")
plt.ylabel("Radius / m")
plt.plot(n[1:-1], r[1:-1])

plt.savefig("Blatt4_Aufgabe3b_Feline,Hannah.pdf")

#print(np.absolute(np.cos(2 * np.pi * 20 * 10**6 * t)))
print(r)
print(g)
print(T)





