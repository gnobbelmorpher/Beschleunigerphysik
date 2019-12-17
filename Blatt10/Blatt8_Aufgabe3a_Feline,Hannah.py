#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, m_e, mu_0
from mpl_toolkits.mplot3d import axes3d
from scipy.stats import norm

#einlesen der txt Datei, in der die Struktur steht

typ, länge, stärke = np.genfromtxt("struktur.txt", skip_header=1,skip_footer=1, unpack=True)
wiederholungen = np.genfromtxt("struktur.txt",skip_footer=typ.size+1, unpack=True)
wiederholungen = int(wiederholungen)

#setzt den zähler bei welchem modul man ist zu beginn auf null (wird später erhöht)
modul=0

#berechnet die zehn matzrizen für den quadrupol

def quadurpol(modul):
    
    #auswahl welche komponente fokussiert wird
    if stärke[modul] < 0: 
        a = 0
        b = 2
    else:
        a = 2
        b = 0

    for i in range(0, 10):


        struktur[a,a, i+modul*10] = np.cos(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a+1,a, i+modul*10] = 1/np.sqrt(np.absolute(stärke[modul])) * np.sin(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a,a+1, i+modul*10] = -np.sqrt(np.absolute(stärke[modul])) * np.sin(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a+1,a+1, i+modul*10] = np.cos(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)

        struktur[b,b, i+modul*10] = np.cosh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b+1,b, i+modul*10] = 1/np.sqrt(np.absolute(stärke[modul])) * np.sinh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b,b+1, i+modul*10] = np.sqrt(np.absolute(stärke[modul])) * np.sinh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b+1,b+1, i+modul*10] = np.cosh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)



#bestimmt die zehn matrizen für die Driftstecke

def driftstrecke(modul):

    
    for i in range(0, 10):


        struktur[1,0, i+modul*10] = länge[modul]/10
        struktur[3,2, i+modul*10] = länge[modul]/10
        


 # n ist modultyp
 # modul ist modulnummer   
 # wählt aus welche art von matrix erzeugt wird, je nach element

def auswahl(n, modul):
    if n == 0:
        driftstrecke(modul)
    if n == 4:
        quadurpol(modul)

# berechnet die teilchenvektoren nach mit den matrizen der Magnetstrukuren

def flug(struktur, teilchen, wiederholungen):
    
    for i in range(0, wiederholungen):
        for n in range(0, 10*(typ.size)):
            teilchen[:,:,n+1+ 10*(typ.size)*i] =  np.matmul( struktur[:,:,n] ,teilchen[:,:,n + 10*(typ.size)*i] )
           

#eigentliches Programm
#prüft zunächst, ob alle einträge gültig sind

if typ[:].all() in [0,1,2,4,6,]:

    #erstellt die matrizen der magnetstrukturen zunächst als einheizmatirzen
    struktur = np.zeros((6,6, 10*(typ.size)))
    for a in range(0, 6):
        struktur[a,a,:]=1

    #überschreiben mit den richtigen matrizen
    for n in typ:
        auswahl(n, modul)
        modul = modul +1

   
    #erstellen der teilchen
    teilchen = np.zeros((6,1000, 10*(typ.size)*wiederholungen+1))
    for a in range(0,1000):
            teilchen[0,a,0] = np.random.normal(0,1e-3)
            teilchen[1,a,0] = np.random.normal(0,1e-3)
            teilchen[2,a,0] = np.random.normal(0,1e-3)
            teilchen[3,a,0] = np.random.normal(0,1e-3)


    flug(struktur, teilchen, wiederholungen)

    #erstellen einer liste mit den orten, an denen matrixeinträge stehen
    s = np.zeros(10*(typ.size)*wiederholungen+1)
    for k in range(0, wiederholungen):
        for m in range(0, typ.size):
            for i in range(0,10):
                s[k*(typ.size)*10 + m*10+i +1] = s[k*(typ.size)*10 + m*10+i] + länge[m]/10


    #plotten
    for i in range(0, 1000):
        plt.plot(s, teilchen[0,i,:])
    
    plt.savefig("Blatt10_Aufgabe3TeilchenspureninX_Feline,Hannah.pdf")

    plt.clf()


    for i in range(0, 1000):
        plt.plot(s, teilchen[2,i,:])
    
    plt.savefig("Blatt10_Aufgabe3TeilchenspureninY_Feline,Hannah.pdf")

    plt.clf()

    for i in range(0, 40):
        plt.subplot(4, 10, i+1)
        plt.plot(teilchen[0,:,i*10], teilchen[1,:,i*10])


    plt.savefig("Blatt10_Aufgabe3PhasenrauminX_Feline,Hannah.pdf")
    
    plt.clf()

    for i in range(0, 40):
        plt.subplot(4, 10, i+1)
        plt.plot(teilchen[2,:,i*10], teilchen[3,:,i*10])


    plt.savefig("Blatt10_Aufgabe3PhasenrauminY_Feline,Hannah.pdf")

    plt.clf()

# durch eine Verringerung der Standadabweichung des Winkels, kommt es vor allem zu einer Verschmalerung der Ellipsen in Phasenraum.
# Zusätzlich Verringern sich die Abweichungen von der Sollbahn ein wenig.
# Entschuldigung, für das kack Format der Plots.


