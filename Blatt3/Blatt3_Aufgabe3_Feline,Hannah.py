#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


#benötigte Funktionen werden definiert

#setzt jedes mal die gitterpunkte der elektroden auf das gewünschte potential, 
#da diese bei der mittelwertbiludung nicht rausgenommen werden
def elektrode():

    gitter[11:18,11:18] = gitter[84:91,84:91] = 1000
    gitter[84:91,11:18] = gitter[11:18,84:91] = -1000

#bildet jeweils die mittelwerte über die nachbarn
#i, j : punkt im gitter an dem wir uns befinden
#x: richtung in die der nachbar in x-richtung ist (-1,+1)
#y: richtung in die der nachbar in y-richtung ist (-1,+1)

#für die ecken 
def mittelwert2(i, j ,x,y):
   return 0.5 * (gitter[i + x, j] + gitter[i, j + y])

#für rand rechts und links
def mittelwert3vert(i, j ,x,):
   return (1/3) * (gitter[i + x, j] + gitter[i, j + 1] + gitter[i, j - 1])

#für rand oben und unten
def mittelwert3hor(i, j ,y,):
   return (1/3) * (gitter[i + 1, j] + gitter[i - 1, j] + gitter[i, j + y])

#für die mitte
def mittelwert4(i, j ,):
   return (1/4) * (gitter[i + 1, j] + gitter[i - 1, j] + gitter[i, j + 1] + gitter[i, j - 1])


#erzeugen der Ausgagssituation
gitter = 2000* np.random.rand(100,100) -1000

elektrode()


for i in range(1,10000):

    #ecken
    gitter[0,0] = mittelwert2(0,0, 1, +1)
    gitter[0,99] = mittelwert2(0,99, 1, -1)
    gitter[99,0] = mittelwert2(99,0, -1, 1)
    gitter[99,99] = mittelwert2(99,99, -1, -1)

    #ränder
    for i in range(1,98):
        gitter[0, i] = mittelwert3vert(0, i , 1)
        gitter[99, i] = mittelwert3vert(99, i , -1)
        gitter[i, 0] = mittelwert3hor(i, 0 , 1)
        gitter[i, 99] = mittelwert3hor(i, 99 , -1)

    #mitte
    for i in range(1,99):
        for j in range(1,99):
            gitter[i,j] = mittelwert4(i, j)
    elektrode()


fig, ax = plt.subplots()
im = ax.imshow(gitter)

cbar = ax.figure.colorbar(im, ax=ax, )
cbar.ax.set_ylabel("elektrisches Potential / V", rotation=-90, va="bottom")

plt.savefig("Blatt3_Aufgabe3_Feline,Hannah.pdf")




