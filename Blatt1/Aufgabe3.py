import numpy as np



#Aufgabeinteil a
text = "Hallo Welt!"

print(text)

file = open("AufgabenteilA.txt", "w")
file.write(text)
file.close


#Aufgabenteil b
zahlen = np.linspace(1, 100, 100)

print(zahlen)

np.savetxt("AufgabenteilB.txt", zahlen, fmt="%d")

#Aufgabenteil c

numbers = np.loadtxt("AufgabenteilB.txt")
prims = []

#für jede zahl in dem eingelesenen array wird eine variable (priemzahl - ja oder nein) auf true gesetzt
#ist die zahl ohne rest durch eine zahl größer 2 und kleiner als sie selbst teilbar wird die variable auf false gesetzt
#ist am ende der schleife die variable noch true, wird die zahl an das array prims angehängt
for number in numbers:
    i = 2
    isPrime = True
    while i < number:
        if number % i == 0:
            isPrime = False

        i += 1
    if isPrime == True:
        prims.append(number)

print(prims)



