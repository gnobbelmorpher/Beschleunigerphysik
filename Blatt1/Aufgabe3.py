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

np.savetxt("AufgabenteilB.txt", zahlen)

#Aufgabenteil c

numbers = np.loadtxt("AufgabenteilB.txt")
prims = []


for number in numbers:
    i = 2
    isPrime = True
    while i < number:
        if number % i == 0:
            isPrime = False

        i += 1
    if isPrime == True:
        prims.append(number)

np.savetxt("AufgabenteilC.txt", prims)



