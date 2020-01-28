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

#von mir ergänzt Blatt 12a
    matrix[0,5]=matrix[4,1]=R*(1-np.cos(L/R))
    matrix[1,5]=matrix[4,0]=np.sin(L/R)
    matrix[4,5]=R*(L/R)-np.sin(L/R)

    return matrix

#Funktion, die die Startwerte der Elektronen würfelt und eine (6 x Anzahl)-große Matrix zurückgibt. Jede Spalte entspricht dem 6er-Vektor eines Elektrons.
#def init_electrons(Anzahl,sigma_x,sigma_x_prime,sigma_y,sigma_y_prime,sigma_z,sigma_E):
#    x=sigma_x*np.random.randn(Anzahl)
#    y=sigma_y*np.random.randn(Anzahl)
#    x_prime=sigma_x_prime*np.random.randn(Anzahl)
#    y_prime=sigma_y_prime*np.random.randn(Anzahl)
#    z=sigma_z*np.random.randn(Anzahl)
#    dEE=sigma_E*np.random.randn(Anzahl)
#    return np.array([x,x_prime,y,y_prime,z,dEE])


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


#Baue alle Matrizen des Rings in ein großes 3D Array. Jede "Scheibe" (alle_Matrizen[:,:,i]) entspricht einer Transfermatrix.
alle_Matrizen=np.zeros([6,6,anzahl_Elemente])
umlauf=np.eye(6)
for i in range(anzahl_Elemente):
    if M_Ring[i]==0:
        alle_Matrizen[:,:,i]=get_Drift(L_Ring[i])
    elif M_Ring[i]==1:
        alle_Matrizen[:,:,i]=get_Edge(L_Ring[i],S_Ring[i])
    elif M_Ring[i]==2:
        alle_Matrizen[:,:,i]=get_Dipole(L_Ring[i],S_Ring[i]) 
    elif M_Ring[i]==4:
        alle_Matrizen[:,:,i]=get_Quad(L_Ring[i],S_Ring[i])
    
        umlauf=np.matmul(alle_Matrizen[:,:,i],umlauf)

#Einträge aus der Umlaufmatrix rausschreiben um optische Funktionen zu berechnen
sxstrich = umlauf[0,0]
sx = - umlauf[0,1]
cxstrich = - umlauf[1,0]
cx = umlauf[1,1]

systrich = umlauf[2,2]
sy = - umlauf[2,3]
cystrich = - umlauf[3,2]
cy = umlauf[3,3]

m13 = umlauf[0,5]
m23 = umlauf[1,5]

betax = np.zeros(anzahl_Elemente+1)
alphax = np.zeros(anzahl_Elemente+1)
gammax = np.zeros(anzahl_Elemente+1)

betay = np.zeros(anzahl_Elemente+1)
alphay = np.zeros(anzahl_Elemente+1)
gammay = np.zeros(anzahl_Elemente+1)

d = np.zeros(anzahl_Elemente+1)
dstrich = np.zeros(anzahl_Elemente+1)
nullarray = np.zeros(anzahl_Elemente+1)

betax[0] = 2*sx/np.sqrt(2-cx**2-2*sx*cxstrich-sxstrich**2)
alphax[0] = (cx - sxstrich)/(2*sx)*betax[0]
gammax[0] = (1-alphax[0]**2)/betax[0]

betay[0] = 2*sy/np.sqrt(2-cy**2-2*sy*cystrich-systrich**2)
alphay[0] = (cy - systrich)/(2*sy)*betay[0]
gammay[0] = (1-alphay[0]**2)/betay[0]

dstrich[0] = (cxstrich*m13+m23*(1-cx))/(2-cx-sxstrich)
d[0]= (sx*dstrich[0]+m13)/(1-cx)

optikx = np.array([betax, alphax, gammax])
optiky = np.array([betay, alphay, gammay])
dispersion = np.array([d, dstrich, nullarray, nullarray, nullarray, np.ones(anzahl_Elemente+1)])
print(dispersion.shape)
print(optikx.shape)

#optische Funktionen nach jedem Element berechnen
for n in range(0, anzahl_Elemente):
    
    transmatrixx= np.array([[alle_Matrizen[1,1,n]**2, 2*alle_Matrizen[0,1,n]*alle_Matrizen[1,1,n], alle_Matrizen[0,1,n]**2],[alle_Matrizen[1,0,n]*alle_Matrizen[1,1,n], alle_Matrizen[0,0,n]*alle_Matrizen[1,1,n]+alle_Matrizen[1,0,n]*alle_Matrizen[0,1,n], alle_Matrizen[0,0,n]*alle_Matrizen[0,1,n]],[alle_Matrizen[1,0,n]**2, 2*alle_Matrizen[0,0,n]*alle_Matrizen[1,0,n], alle_Matrizen[0,0,n]**2]])
    transmatrixy= np.array([[alle_Matrizen[3,3,n]**2, 2*alle_Matrizen[2,3,n]*alle_Matrizen[3,3,n], alle_Matrizen[2,3,n]**2],[alle_Matrizen[3,2,n]*alle_Matrizen[3,3,n], alle_Matrizen[2,2,n]*alle_Matrizen[3,3,n]+alle_Matrizen[3,2,n]*alle_Matrizen[2,3,n], alle_Matrizen[2,2,n]*alle_Matrizen[2,3,n]],[alle_Matrizen[3,2,n]**2, 2*alle_Matrizen[2,2,n]*alle_Matrizen[3,2,n], alle_Matrizen[2,2,n]**2]])

    optikx[:,n+1]=np.dot(transmatrixx[:,:],optikx[:,n])
    optiky[:,n+1]=np.dot(transmatrixy[:,:],optiky[:,n])
    dispersion[:, n+1] = np.dot(alle_Matrizen[:,:,n],dispersion[:, n])

#phasenvorschub phix und momentum compaction factor a
phix = 0
for beta in optikx[0,:]:
    phix = phix + 1/beta

print(phix)

phiy = 0
for beta in optiky[0,:]:
    phiy = phiy + 1/beta

print(phiy)

a = 0
for d in dispersion[0,:]:
    a = d/3.82
a = a / s[-1]
print(a)

#Elektronen generieren und in die "vorderste Scheibe" eines 3D Arrays packen. 
#N_elec=1000
#sigma_x=1e-3
#sigma_x_Strich=1e-3
#sigma_y=1e-3
#sigma_y_Strich=1e-3
#sigma_z=0
#sigma_dEE=0

#elecs=np.zeros([6,N_elec,anzahl_Elemente+1])    #3D Array erstellen, mit anzahl_Elemente+1 "Scheiben" (Startwerte + Werte nach jedem Element)
#elecs[:,:,0]=init_electrons(N_elec,sigma_x,sigma_x_Strich,sigma_y,sigma_y_Strich,sigma_z,sigma_dEE)

#Schleife über alle Elemente um 
#for i in range(anzahl_Elemente):
#    elecs[:,:,i+1]=np.dot(alle_Matrizen[:,:,i],elecs[:,:,i])
    
############    
## Plotten##
############
    
#Plot der x und y Ablagen gegen s
#fig1, (fig1_ax1,fig1_ax2)=plt.subplots(2,1)
#fig1_ax1.plot(s,elecs[0,:,:].T,'-',lw=0.1,c='tab:blue')
#fig1_ax1.set_xlabel('longtidinale Postion s in m')
#fig1_ax1.set_ylabel('horizontale Ablage x in m')
#fig1_ax1.set_ylim([-0.025,0.025])
#fig1_ax2.plot(s,elecs[2,:,:].T,'-',lw=0.1,c='tab:blue')
#fig1_ax2.set_xlabel('longtidinale Postion s in m')
#fig1_ax2.set_ylabel('vertikale Ablage y in m')
#fig1_ax2.set_ylim([-0.025,0.025])

#plt.subplots_adjust(hspace=0.4)

#Plot der optischen Funktionen gegen s
fig3, (fig3_ax1,fig3_ax2, fig3_ax3, fig3_ax4, fig3_ax5)=plt.subplots(5,1)
fig3_ax1.plot(s,optikx[0,:].T,'-',c='tab:blue')
fig3_ax1.set_xlabel('longtidinale Postion s in m')
fig3_ax1.set_ylabel('Beta-Funktion in x')
#fig3_ax1.set_ylim([-0.025,0.025])
fig3_ax2.plot(s,optikx[1,:].T,'-',c='tab:blue')
fig3_ax2.set_xlabel('longtidinale Postion s in m')
fig3_ax2.set_ylabel('Alpha-Funktion')
#fig3_ax2.set_ylim([-0.025,0.025])
fig3_ax3.plot(s,optikx[2,:].T,'-',c='tab:blue')
fig3_ax3.set_xlabel('longtidinale Postion s in m')
fig3_ax3.set_ylabel('Gamma-Funktion')
#ig3_ax3.set_ylim([-0.025,0.025])
fig3_ax4.plot(s,dispersion[0,:].T,'-',c='tab:blue')
fig3_ax4.set_xlabel('longtidinale Postion s in m')
fig3_ax4.set_ylabel('Dispersion')
#ig3_ax4.set_ylim([-0.025,0.025])
fig3_ax5.plot(s,dispersion[1,:].T,'-',c='tab:blue')
fig3_ax5.set_xlabel('longtidinale Postion s in m')
fig3_ax5.set_ylabel('Ableitung der Dispersion')
#ig3_ax5.set_ylim([-0.025,0.025])

plt.subplots_adjust(hspace=0.4)

## Animation des Phasenraums entlang des Rings
#fig2, (fig2_ax1,fig2_ax2)=plt.subplots(1,2)

#fig2_ax1.set_aspect(1) 
#fig2_ax1.set_xlim([-0.025,0.025])
#fig2_ax1.set_ylim([-0.025,0.025])
#fig2_ax1.set_xlabel('horizontale Ablabe x in mm')
#fig2_ax1.set_ylabel("horizontaler Winkel x' in mrad")
#phasespaceX, = fig2_ax1.plot([],[],'.',ms=0.5)

#fig2_ax2.set_aspect(1) 
#fig2_ax2.set_xlim([-0.025,0.025])
#fig2_ax2.set_ylim([-0.025,0.025])
#fig2_ax2.set_xlabel('vertikal Ablabe y in mm')
#fig2_ax2.set_ylabel("vertikaler Winkel y' in mrad")
#phasespaceY, = fig2_ax2.plot([],[],'.',ms=0.5)
#plt.subplots_adjust(wspace=0.5)
#line=[phasespaceX,phasespaceY]

#def init_animation():
#    line[0].set_data([],[])
#    line[1].set_data([],[])
#    return line    

#def animate(i):
#    line[0].set_data(elecs[0,:,i],elecs[1,:,i])
#    line[1].set_data(elecs[2,:,i],elecs[3,:,i])
#    if M_Ring[i] == 0:
#        Element = 'Driftstrecke'
#    if M_Ring[i] == 1:
#        Element = 'Kante'
#    if M_Ring[i] == 2:
#        Element = 'schwache Fokussierung'
#    if M_Ring[i] == 4:
#        Element = 'Quadrupol'
#    fig2_ax1.set_title('s = %6.2f\n' % s[i] + Element)
#    fig2_ax2.set_title('s = %6.2f\n' % s[i] + Element)
#    return line

#anim=animation.FuncAnimation(fig2,animate,range(anzahl_Elemente), init_func=init_animation, interval=100)
plt.show()