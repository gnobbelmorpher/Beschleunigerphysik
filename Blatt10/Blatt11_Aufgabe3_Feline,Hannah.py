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
        struktur[a,a+1, i+modul*10] = 1/np.sqrt(np.absolute(stärke[modul])) * np.sin(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a+1,a, i+modul*10] = -np.sqrt(np.absolute(stärke[modul])) * np.sin(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[a+1,a+1, i+modul*10] = np.cos(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)

        struktur[b,b, i+modul*10] = np.cosh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b,b+1, i+modul*10] = 1/np.sqrt(np.absolute(stärke[modul])) * np.sinh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b+1,b, i+modul*10] = np.sqrt(np.absolute(stärke[modul])) * np.sinh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)
        struktur[b+1,b+1, i+modul*10] = np.cosh(np.sqrt(np.absolute(stärke[modul]))*länge[modul]/10)



#bestimmt die zehn matrizen für die Driftstecke

def driftstrecke(modul):

    
    for i in range(0, 10):


        struktur[0,1, i+modul*10] = länge[modul]/10
        struktur[2,3, i+modul*10] = länge[modul]/10
        


 # n ist modultyp
 # modul ist modulnummer   
 # wählt aus welche art von matrix erzeugt wird, je nach element

def dipol(modul):

    #R=U/2pi=länge*anzahl/2pi=1.5*16/2pi=3,8197=stärke

        for i in range(0, 10):

            struktur[0,0, i+modul*10] = np.cos(1/stärke[modul]*länge[modul]/10)
            struktur[0,1, i+modul*10] = stärke[modul] * np.sin(1/stärke[modul]*länge[modul]/10)
            struktur[1,0, i+modul*10] = -1/stärke[modul] * np.sin(1/stärke[modul]*länge[modul]/10)
            struktur[1,1, i+modul*10] = np.cos(1/stärke[modul]*länge[modul]/10)

            struktur[2,3, i+modul*10] = länge[modul]/10

            struktur[0,5, i+modul*10] = stärke[modul] * (1- np.cos(länge[modul]/stärke[modul]))
            struktur[1,5, i+modul*10] = np.sin(länge[modul]/stärke[modul])
            struktur[4,0, i+modul*10] = np.sin(länge[modul]/stärke[modul])
            struktur[4,1, i+modul*10] = stärke[modul] * (1- np.cos(länge[modul]/stärke[modul]))
            struktur[4,5, i+modul*10] = stärke[modul] * (länge[modul]/stärke[modul] - np.sin(länge[modul]/stärke[modul]))


def kante(modul):

    if typ[modul-1] == 2:
        struktur[1,0, modul*10] = np.tan(länge[modul-1]/2*stärke[modul-1])/stärke[modul-1]
        struktur[4,3, modul*10] = -np.tan(länge[modul-1]/2*stärke[modul-1])/stärke[modul-1]
    else:
        struktur[1,0, modul*10] = np.tan(länge[modul+1]/2*stärke[modul+1])/stärke[modul+1]
        struktur[4,3, modul*10] = -np.tan(länge[modul+1]/2*stärke[modul+1])/stärke[modul+1]

def auswahl(n, modul):
    if n == 0:
        driftstrecke(modul)
    if n == 4:
        quadurpol(modul)
    if n == 2:
        dipol(modul)
    if n == 1:
        kante(modul)

# berechnet die teilchenvektoren nach mit den matrizen der Magnetstrukuren

def flug(struktur, teilchen, wiederholungen):

    # Blatt 12
# flexibler machen für andere Strukturen!

    for i in range(0, 8):
        for k in range(0, 129):
            umlauf = np.matmul(struktur[:,:,k +1], struktur[:,:,k])
    print(umlauf)

    sxstrich = umlauf[0,0]
    sx = - umlauf[0,1]
    cxstrich = - umlauf[1,0]
    cx = umlauf[1,1]

    systrich = umlauf[2,2]
    sy = - umlauf[2,3]
    cystrich = - umlauf[3,2]
    cy = umlauf[3,3]

    betax = np.zeros(130)
    alphax = np.zeros(130)
    gammax = np.zeros(130)

    betay = np.zeros(130)
    alphay = np.zeros(130)
    gammay = np.zeros(130)

    betax[0] = 2*sx/np.sqrt(2-cx**2-2*sx*cxstrich-sxstrich**2)
    alphax[0] = (cx - sxstrich)/(2*sx)*betax[0]
    gammax[0] = (1-alphax[0]**2)/betax[0]

    betay[0] = 2*sy/np.sqrt(2-cy**2-2*sy*cystrich-systrich**2)
    alphay[0] = (cy - systrich)/(2*sy)*betay[0]
    gammay[0] = (1-alphay[0]**2)/betay[0]

    optikx = np.array([betax, alphax, gammax])
    optiky = np.array([betay, alphay, gammay])
    print(optikx)
    # irgendwo muss ein Vorzeichenfehler sein, weil eine Wurzel negativ wird
    #Bei der Dispersion bin ich verwirrt welche Matrixeinträge jetzt relevant sind bzw welche Matrix

    
    
    for i in range(0, wiederholungen):
        for n in range(0, 10*(typ.size)):
            teilchen[:,:,n+1+ 10*(typ.size)*i] =  np.matmul( struktur[:,:,n] ,teilchen[:,:,n + 10*(typ.size)*i] )

            transmatrixx= [[struktur[1,1,n]**2, 2*struktur[0,1,n]*struktur[1,1,n], struktur[0,1,n]**2],
            [struktur[1,0,n]*struktur[1,1,n], struktur[0,0,n]*struktur[1,1,n]+struktur[1,0,n]*struktur[0,1,n], struktur[0,0,n]*struktur[0,1,n]],
            [struktur[1,0,n]**2, 2*struktur[0,0,n]*struktur[1,0,n], struktur[0,0,n]**2]]
            transmatrixy= [[struktur[3,3,n]**2, 2*struktur[2,3,n]*struktur[3,3,n], struktur[2,3,n]**2],
            [struktur[3,2,n]*struktur[3,3,n], struktur[2,2,n]*struktur[3,3,n]+struktur[3,2,n]*struktur[2,3,n], struktur[2,2,n]*struktur[2,3,n]],
            [struktur[3,2,n]**2, 2*struktur[2,2,n]*struktur[3,2,n], struktur[2,2,n]**2]]

            optikx[n+1+ 10*(typ.size)*i,:] =  np.matmul( transmatrixx[:,:,n] ,optikx[n + 10*(typ.size)*i,:] )
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
    
    plt.savefig("Blatt11_Aufgabe3TeilchenspureninX_Feline,Hannah.pdf")

    plt.clf()


    for i in range(0, 1000):
        plt.plot(s, teilchen[2,i,:])
    
    plt.savefig("Blatt11_Aufgabe3TeilchenspureninY_Feline,Hannah.pdf")

    plt.clf()

    for i in range(0, 40):
        plt.subplot(4, 10, i+1)
        plt.plot(teilchen[0,:,i*10], teilchen[1,:,i*10])


    plt.savefig("Blatt11_Aufgabe3PhasenrauminX_Feline,Hannah.pdf")
    
    plt.clf()

    for i in range(0, 40):
        plt.subplot(4, 10, i+1)
        plt.plot(teilchen[2,:,i*10], teilchen[3,:,i*10])


    plt.savefig("Blatt11_Aufgabe3PhasenrauminY_Feline,Hannah.pdf")

    plt.clf()

# durch eine Verringerung der Standadabweichung des Winkels, kommt es vor allem zu einer Verschmalerung der Ellipsen in Phasenraum.
# Zusätzlich Verringern sich die Abweichungen von der Sollbahn ein wenig.
# Entschuldigung, für das kack Format der Plots.


