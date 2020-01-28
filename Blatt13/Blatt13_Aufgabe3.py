import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

#########################
##Funktionen definieren##
#########################

# Funktion, um das Lattice-File einzulesen und einige fehlerhafte Eingaben abzufangen.
def read_lattice(filename):
    with open(filename,'r') as f:
        lines=f.readlines()
        
    if not lines[0][:-1].isdigit():
        raise ValueError('Das File beginnt nicht mit einem Integer für N!')
        
    if lines[-1]!='-1':
        raise ValueError('Das Latticefile endet nicht mit -1! (Ggf. überprüfen, ob am Ende eine Leerzeile steht)')
        
    N=int(lines[0][:-1])
    M=np.zeros(np.size(lines)-2)
    L=np.zeros(np.size(lines)-2)
    S=np.zeros(np.size(lines)-2)
    k=0
    for line in lines[1:-1]:
        M[k],L[k]=line.split()[0:2]
        if M[k] not in [0,1,2,4,6]:
            raise ValueError('Ringelement nicht erkannt!')
        if M[k] != 0:
            S[k]=line.split()[2]
        k+=1
    return M,L,S,N

# Funktion, die eine Transfermatrix für eine Driftstrecke der Länge L zurück gibt.
def get_Drift(L):
    matrix=np.eye(6)
    matrix[0,1]=matrix[2,3]=L
    return matrix

# Funktion, die eine Transfermatrix für einen Quadrupol der Länge L und Stärke k zurück gibt.
def get_Quad(L,k):
    omega=np.sqrt(np.abs(k))*L
    defokussierend = np.array([[np.cosh(omega),1/np.sqrt(np.abs(k))*np.sinh(omega)],[np.sqrt(np.abs(k))*np.sinh(omega),np.cosh(omega)]])
    fokussierend = np.array([[np.cos(omega),1/np.sqrt(np.abs(k))*np.sin(omega)],[-np.sqrt(np.abs(k))*np.sin(omega),np.cos(omega)]])
    if k>0:
        matrix=np.identity(6)
        matrix[0:2,0:2]=defokussierend
        matrix[2:4,2:4]=fokussierend
    elif k<0:
        matrix=np.identity(6)
        matrix[0:2,0:2]=fokussierend
        matrix[2:4,2:4]=defokussierend
    else:
        matrix=get_Drift(L)  #falls ein Quadrupol ausgeschaltet sein soll (k=0), geht die Quadrupolmatrix in eine Driftmatrix über.
    return matrix


# Funktion, die eine Transfermatrix für eine Kante eines Dipols mit Länge L und Biegeradius R zurück gibt.
def get_Edge(L,R):
    Psi=L/(2*R) #halber Biegewinkel
    matrix=np.eye(6)
    matrix[1,0]=np.tan(Psi)/R
    matrix[3,2]=-np.tan(Psi)/R
    return matrix

# Funktion, die eine Transfermatrix eines Dipols mit Länge L und Biegeradius R zurück gibt.
def get_Dipole(L,R):
    Psi=L/R     #Biegewinkel
    matrix=np.eye(6)
    matrix[0,0]=matrix[1,1]=np.cos(Psi)
    matrix[0,1]=R*np.sin(Psi)
    matrix[1,0]=-np.sin(Psi)/R
    matrix[2,3]=L
    matrix[0,5]=R*(1-np.cos(L/R))
    matrix[1,5]=np.sin(L/R)
    matrix[4,0]=np.sin(L/R)
    matrix[4,1]=R*(1-np.cos(L/R))
    matrix[4,5]=R*(L/R-np.sin(L/R))
    return matrix

#Funktion, die die Startwerte der Elektronen würfelt und eine (6 x Anzahl)-große Matrix zurückgibt. Jede Spalte entspricht dem 6er-Vektor eines Elektrons.
def init_electrons(Anzahl,sigma_x,sigma_x_prime,sigma_y,sigma_y_prime,sigma_z,sigma_E):
    x=sigma_x*np.random.randn(Anzahl)
    y=sigma_y*np.random.randn(Anzahl)
    x_prime=sigma_x_prime*np.random.randn(Anzahl)
    y_prime=sigma_y_prime*np.random.randn(Anzahl)
    z=sigma_z*np.random.randn(Anzahl)
    dEE=sigma_E*np.random.randn(Anzahl)
    return np.array([x,x_prime,y,y_prime,z,dEE])

def get_OneTurn(matrizen,N):
    oneTurn=np.eye(6)
    for i in range(N):
            oneTurn=np.dot(matrizen[:,:,i],oneTurn)
    return oneTurn

def get_alleTwiss(matrizen,N):
    alleTwiss_x=np.zeros([3,3,N])
    alleTwiss_y=np.zeros([3,3,N])
    for i in range(N):
        matrix=matrizen[:,:,i]
        C_x=matrix[1,1]
        S_x=-matrix[0,1]
        C_prime_x=-matrix[1,0]
        S_prime_x=matrix[0,0]
        
        C_y=matrix[3,3]
        S_y=-matrix[2,3]
        C_prime_y=-matrix[3,2]
        S_prime_y=matrix[2,2]
        
        alleTwiss_x[:,:,i]=np.array([[C_x**2,-2*S_x*C_x,S_x**2],[-C_x*C_prime_x,S_x*C_prime_x+S_prime_x*C_x,-S_x*S_prime_x],[C_prime_x**2,-2*S_prime_x*C_prime_x,S_prime_x**2]])
        alleTwiss_y[:,:,i]=np.array([[C_y**2,-2*S_y*C_y,S_y**2],[-C_y*C_prime_y,S_y*C_prime_y+S_prime_y*C_y,-S_y*S_prime_y],[C_prime_y**2,-2*S_prime_y*C_prime_y,S_prime_y**2]])
    return alleTwiss_x,alleTwiss_y
#######################################
## Beginn der eigentlichen Berechnung##
#######################################
    

#Magnetstruktur einlesen
M,L,S,N=read_lattice('lattice.txt')

#Neue M,L,S Vektoren bauen, in denen die Daten des gesamten Rings stehen, inklusive der Teilung der Magnete.
#Außerdem, den Ort jeden Elements (genau genommen das Ende eines jeden Elements) in den Vektor s schreiben.
Teilung=10      #In wieviele Teile ein Magnet gestückelt wird
M_Ring=[]
L_Ring=[]
S_Ring=[]
s=[0]       #beginnt mit 0, der Startposition aller Elektronen

for k in range(N):                          #Folgendes wird N-mal ausgeführt
    for i in range(len(M)):                 #Schleife über die Elemente im File
        if M[i]==1:                         #Kante abfangen, da diese keine Länge hat brauch sie auch nicht geteilt werden
            M_Ring.append(M[i])
            L_Ring.append(L[i])
            S_Ring.append(S[i])
            s.append(s[-1])                 #Kante hat keine Länge, daher hier vorherige Position beibehalten
            continue                        #Überspringen der Teilungsschleife
        else:            
            for j in range(Teilung):            #sorgt für die Kleinteilung der Magnete
                M_Ring.append(M[i])
                L_Ring.append(L[i]/Teilung)
                S_Ring.append(S[i])
                s.append(s[-1]+L[i]/Teilung)

anzahl_Elemente=len(M_Ring) #Gesamtzahl aller Elemente im Ring (inklusive der Teilung der Magnete und der Periodizität gegeben mit N)

#Parameter zur Berechnung der Chromatizität
C=np.zeros(anzahl_Elemente+1)
#Baue alle Matrizen des Rings in ein großes 3D Array. Jede "Scheibe" (alle_Matrizen[:,:,i]) entspricht einer Transfermatrix.
alle_Matrizen=np.zeros([6,6,anzahl_Elemente])
for i in range(anzahl_Elemente):
    if M_Ring[i]==0:
        alle_Matrizen[:,:,i]=get_Drift(L_Ring[i])
    elif M_Ring[i]==1:
        alle_Matrizen[:,:,i]=get_Edge(L_Ring[i],S_Ring[i])
    elif M_Ring[i]==2:
        alle_Matrizen[:,:,i]=get_Dipole(L_Ring[i],S_Ring[i]) 
        C[i]=S_Ring[i]
    elif M_Ring[i]==4:
        alle_Matrizen[:,:,i]=get_Quad(L_Ring[i],S_Ring[i])
        C[i]=S_Ring[i]
    elif M_Ring[i]==6:
        alle_Matrizen[:,:,i]=get_Drift(L_Ring[i])

# alle Transfermatrizen für die optischen Funktionen generieren
alleTwiss_x,alleTwiss_y=get_alleTwiss(alle_Matrizen,anzahl_Elemente)

#OneTurnMatrix erstellen
OneTurnMatrix=get_OneTurn(alle_Matrizen,anzahl_Elemente)

#Dispersion bei s=0 berechnen
D_prime_0=(OneTurnMatrix[1,0]*OneTurnMatrix[0,5]+OneTurnMatrix[1,5]*(1-OneTurnMatrix[0,0]))/(2-OneTurnMatrix[0,0]-OneTurnMatrix[1,1])
D_0=(OneTurnMatrix[0,1]*D_prime_0+OneTurnMatrix[0,5])/(1-OneTurnMatrix[0,0])

#Startwerte der optischen Funktionen aus der OneTurnMatrix ziehen und in einen Vektor schreiben

C_x=OneTurnMatrix[0,0]
S_x=OneTurnMatrix[0,1]
C_prime_x=OneTurnMatrix[1,0]
S_prime_x=OneTurnMatrix[1,1]

C_y=OneTurnMatrix[2,2]
S_y=OneTurnMatrix[2,3]
C_prime_y=OneTurnMatrix[3,2]
S_prime_y=OneTurnMatrix[3,3]

beta_x_0=np.abs(2*S_x/np.sqrt(2-C_x**2-2*S_x*C_prime_x-S_prime_x**2))
alpha_x_0=(C_x-S_prime_x)*beta_x_0/(2*S_x)
gamma_x_0=(1+alpha_x_0**2)/beta_x_0

beta_y_0=np.abs(2*S_y/np.sqrt(2-C_y**2-2*S_y*C_prime_y-S_prime_y**2))
alpha_y_0=(C_y-S_prime_y)*beta_y_0/(2*S_y)
gamma_y_0=(1+alpha_y_0**2)/beta_y_0

twiss0_x=np.array([beta_x_0,alpha_x_0,gamma_x_0])
twiss0_y=np.array([beta_y_0,alpha_y_0,gamma_y_0])

#Dipersion und optische Funktionen an jedem Element berechnen
D_vec=np.zeros([6,anzahl_Elemente+1])
D_vec[:,0]=np.array([D_0,D_prime_0,0,0,0,1])

twiss_x=np.zeros([3,anzahl_Elemente+1])
twiss_y=np.zeros([3,anzahl_Elemente+1])

twiss_x[:,0]=twiss0_x
twiss_y[:,0]=twiss0_y

for i in range(anzahl_Elemente):
    twiss_x[:,i+1]=np.matmul(alleTwiss_x[:,:,i],twiss_x[:,i])
    twiss_y[:,i+1]=np.matmul(alleTwiss_y[:,:,i],twiss_y[:,i])
    D_vec[:,i+1]=np.matmul(alle_Matrizen[:,:,i],D_vec[:,i])

#momentum compaction factor:
Cy=np.copy(C)
D=D_vec[0,:]
R=np.ones(anzahl_Elemente+1)*1e50
for i in range(anzahl_Elemente):
    if M_Ring[i]==2:
        R[i]=S_Ring[i]
        C[i]=-C[i]*D_vec[0,i]
        Cy[i]=0

alpha=1/s[-1] * np.trapz(D/R,s)   

#Phasenvorschub
QX=1/(2*np.pi)*np.trapz(1/twiss_x[0,:],s)
QY=1/(2*np.pi)*np.trapz(1/twiss_y[0,:],s)

#Chromazität
xiX=1/(4*np.pi)*np.trapz(twiss_x[0,:]*C,s)
xiY=-1/(4*np.pi)*np.trapz(twiss_y[0,:]*Cy,s)

#Emittanz
lorenz=1.7e9 /511e3 +1
c_lorenz= 3.82e-13
a = 0
for i in range(anzahl_Elemente):
        if M_Ring[i]==2:
            epsilon= c_lorenz*lorenz**2/(S_Ring[i]*10*L_Ring[i])*np.trapz(twiss_x[2,i:i+10]*D_vec[0,i:i+10]**2+2*twiss_x[1,i:i+10]*D_vec[1,i:i+10]*D_vec[0,i:i+10]+twiss_x[0,i:i+10]*D_vec[1,i:i+10]**2,s[i:i+10])
            break



print('momentum compaction factor: '+str(alpha))
print('hor. Arbeitspunkt: '+str(QX))
print('ver. Arbeitspunkt: '+str(QY))
print('hor. Chromatizität: '+str(xiX))
print('ver. Chromatizität: '+str(xiY))
print('hor. Emittanz: '+str(epsilon))

plt.figure()
plt.plot(s,twiss_x[0,:],label=r'$\beta_x$')
plt.plot(s,twiss_x[1,:],label=r'$\alpha_x$')
plt.plot(s,twiss_x[2,:],label=r'$\gamma_x$')

plt.plot(s,twiss_y[0,:],label=r'$\beta_y$')
plt.plot(s,twiss_y[1,:],label=r'$\alpha_y$')
plt.plot(s,twiss_y[2,:],label=r'$\gamma_y$')
plt.legend()

plt.figure()
plt.plot(s,D_vec[0,:],label=r'Dispersion $D$')
plt.plot(s,D_vec[1,:],label='D\'')
plt.legend()

#Elektronen generieren und in die "vorderste Scheibe" eines 3D Arrays packen. 
N_elec=1000
sigma_x=1e-3
sigma_x_Strich=1e-3
sigma_y=1e-3
sigma_y_Strich=1e-3
sigma_z=0
sigma_dEE=0

elecs=np.zeros([6,N_elec,anzahl_Elemente+1])    #3D Array erstellen, mit anzahl_Elemente+1 "Scheiben" (Startwerte + Werte nach jedem Element)
elecs[:,:,0]=init_electrons(N_elec,sigma_x,sigma_x_Strich,sigma_y,sigma_y_Strich,sigma_z,sigma_dEE)

#Schleife über alle Elemente um 
for i in range(anzahl_Elemente):
    elecs[:,:,i+1]=np.dot(alle_Matrizen[:,:,i],elecs[:,:,i])
    
#############    
### Plotten##
#############
    
#Plot der x und y Ablagen gegen s
fig1, (fig1_ax1,fig1_ax2)=plt.subplots(2,1)
fig1_ax1.plot(s,elecs[0,:,:].T,'-',lw=0.1,c='tab:blue')
fig1_ax1.set_xlabel('longtidinale Postion s in m')
fig1_ax1.set_ylabel('horizontale Ablage x in m')
fig1_ax1.set_ylim([-0.025,0.025])
fig1_ax2.plot(s,elecs[2,:,:].T,'-',lw=0.1,c='tab:blue')
fig1_ax2.set_xlabel('longtidinale Postion s in m')
fig1_ax2.set_ylabel('vertikale Ablage y in m')
fig1_ax2.set_ylim([-0.025,0.025])

plt.subplots_adjust(hspace=0.4)
## Animation des Phasenraums entlang des Rings
fig2, (fig2_ax1,fig2_ax2)=plt.subplots(1,2)

fig2_ax1.set_aspect(1) 
fig2_ax1.set_xlim([-0.025,0.025])
fig2_ax1.set_ylim([-0.025,0.025])
fig2_ax1.set_xlabel('horizontale Ablabe x in mm')
fig2_ax1.set_ylabel("horizontaler Winkel x' in mrad")
phasespaceX, = fig2_ax1.plot([],[],'.',ms=0.5)

fig2_ax2.set_aspect(1) 
fig2_ax2.set_xlim([-0.025,0.025])
fig2_ax2.set_ylim([-0.025,0.025])
fig2_ax2.set_xlabel('vertikal Ablabe y in mm')
fig2_ax2.set_ylabel("vertikaler Winkel y' in mrad")
phasespaceY, = fig2_ax2.plot([],[],'.',ms=0.5)
plt.subplots_adjust(wspace=0.5)
line=[phasespaceX,phasespaceY]

def init_animation():
    line[0].set_data([],[])
    line[1].set_data([],[])
    return line    

def animate(i):
    line[0].set_data(elecs[0,:,i],elecs[1,:,i])
    line[1].set_data(elecs[2,:,i],elecs[3,:,i])
    if M_Ring[i] == 0:
        Element = 'Driftstrecke'
    if M_Ring[i] == 1:
        Element = 'Kante'
    if M_Ring[i] == 2:
        Element = 'schwache Fokussierung'
    if M_Ring[i] == 4:
        Element = 'Quadrupol'
    if M_Ring[i] == 6:
        Element = 'Sextupol'
    fig2_ax1.set_title('s = %6.2f\n' % s[i] + Element)
    fig2_ax2.set_title('s = %6.2f\n' % s[i] + Element)
    return line

anim=animation.FuncAnimation(fig2,animate,range(anzahl_Elemente), init_func=init_animation, interval=100)
#plt.show()