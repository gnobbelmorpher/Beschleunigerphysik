#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e, mu_0
from mpl_toolkits.mplot3d import axes3d
from scipy.stats import norm

modul=0

typ, länge, stärke = np.genfromtxt("struktur.txt", skip_header=1,skip_footer=1, unpack=True)
wiederholungen = np.genfromtxt("struktur.txt",skip_footer=typ.size+1, unpack=True)
wiederholungen = int(wiederholungen)
print(wiederholungen)

def quadurpol(modul):
    
    if stärke[modul] < 0: 
        a = 0
        b = 2
    else:
        a = 2
        b = 0

    for i in range(0, 10):


        struktur[a,a, i+modul*10] = np.cos(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a+1,a, i+modul*10] = 1/np.sqrt(np.absolute(stärke[modul])) * np.sin(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a,a+1, i+modul*10] = np.sqrt(np.absolute(stärke[modul])) * np.sin(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a+1,a+1, i+modul*10] = np.cos(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)

        struktur[b,b, i+modul*10] = np.cosh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b+1,b, i+modul*10] = 1/np.sqrt(np.absolute(stärke[modul])) * np.sinh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b,b+1, i+modul*10] = np.sqrt(np.absolute(stärke[modul])) * np.sinh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b+1,b+1, i+modul*10] = np.cosh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)

def driftstrecke(modul):

    
    for i in range(0, 10):


        struktur[1,0, i+modul*10] = länge[modul]/10
        struktur[3,2, i+modul*10] = länge[modul]/10
        
 # n ist modultyp
 # modul ist modulnummer   

def auswahl(n, modul):
    if n == 0:
        driftstrecke(modul)
    if n == 4:
        quadurpol(modul)

def flug(struktur, teilchen, wiederholungen):
    
    for i in range(0, wiederholungen):
        print(i)
        for n in range(0, 10):
           # teilchen[:,:,n+1+ 10*(typ.size)*i] =  np.matmul( struktur[:,:,n] ,teilchen[:,:,n + 10*(typ.size)*i] )
           print(teilchen[:, :, :])
           


if typ[:].all() in [0,1,2,4,6,]:
    struktur = np.zeros((6,6, 10*(typ.size)))
    for a in range(0, 6):
        struktur[a,a,:]=1
    for n in typ:
        auswahl(n, modul)
        modul = modul +1
    
    teilchen = np.zeros((1000,6, 10*(typ.size)*wiederholungen))
    for a in range(0,1000):
            teilchen[a,0,0] = np.random.normal(0,1)
            teilchen[a,1,0] = np.random.normal(0,1)
            teilchen[a,2,0] = np.random.normal(0,1)
            teilchen[a,3,0] = np.random.normal(0,1)

    flug(struktur, modul, wiederholungen)

#plt.title("Teilchenbewegung im longitudinalen Phasenraum")

#plt.xlabel(r"$\Delta \psi / rad$")
#plt.ylabel(r"$\Delta E / E$")

#plt.savefig("Blatt9_Aufgabe3a_Feline,Hannah.pdf")


#stabile Trajektorien: Phase: ca. -2,5 bis 2,5
#                      Energie: bis Abweichung von 0,34 stabil
# Elektronen im instabilen Bereich verlieren immer weiter an Energie 



