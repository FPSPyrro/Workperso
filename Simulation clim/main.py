import math
from math import *
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import *

#-------------------
#Paramètres fixes
#-------------------
FichierMaillage="circle" #nom du fichier (pas besoin de placé .msh juste le nom est nécessaire)
Kelvin = 273.15
c = 4.2 * 10**6 #CONSTANTE
Keau = 0.6 #w M-1 k-1
f = 0
R = 0.05 # m
L = 50 # m
Ue = 40+Kelvin # °c
n = 10 #itération dans le schéma

#-------------------
#Paramètres de réglages
#-------------------
V = 1.4/3.6 # 0.5 km/h soit V = 0.5/3.6 #m/s
alpha = 26/0.002 # Kinox / Einox
Uc = 14+Kelvin #°c  

#-------------------
#Paramètres iter-dépendants
#-------------------
T = L/V #temps de parcours du flot
delta = T/n

def lit_fichier_msh(fileName):
    f = open("msh/"+fileName+".msh", "r")
    #nbn,nbe,nba = int
    #coord,tri,ar,refn,reft,refa = [],[],[],[],[],[],[],[],[]

    line = f.readline()
    data = line.split()

    nbn = int(data[0])
    nbe = int(data[1])
    nba = int(data[2])
    coord = zeros((nbn,2))
    #coord = [[0 for i in range(2)] for j in range(nbn)]
    #refn = [0 for j in range(nbn)]
    refn = zeros(nbn)
    for i in range(nbn):
        line = f.readline()
        data = line.split()

        
        coord[i][0] = double(data[0])
        coord[i][1] = double(data[1])
        refn[i] = int(data[2])
    #print(coord,refn)
    tri = zeros((nbe,3))
    #tri = [[0 for i in range(3)] for j in range(nbe)]
    #reft = [0 for j in range(nbe)]
    reft = zeros(nbe)
    for i in range(nbe):
        line = f.readline()
        data = line.split()

        
        tri[i][0] = int(data[0])-1
        tri[i][1] = int(data[1])-1
        tri[i][2] = int(data[2])-1
        reft[i] = data[3]
    #print(tri,reft,nbe)
    ar = zeros((nba,2))
    #ar = [[0 for i in range(2)] for j in range(nba)]
    #refa = [0 for j in range(nba)]
    refa = zeros(nba)
    for i in range(nba):
        line = f.readline()
        data = line.split()

        
        ar[i][0] = int(data[0])-1
        ar[i][1] = int(data[1])-1
        refa[i] = data[2]
    #print(ar,refa)
    f.close()
    return nbn,nbe,nba,coord,tri,ar,refn,reft,refa
def pas_et_qualite_triangle(s1,s2,s3):
    la = zeros(3)
    la[0] = linalg.norm(s2-s1)
    la[1] = linalg.norm(s3-s1)
    la[2] = linalg.norm(s2-s3)
    pas = max(la)
    dt = 0.5*sum(la)
    aire = 0.5*abs(cross(s2-s1,s3-s1))
    rayon = aire/dt
    qualite = (sqrt(3)/6) * (pas/rayon)
    return pas,qualite
def pas_et_qualite_maillage(nbe,coord,tri):
    maxPas = 0.0
    maxQual = 1.0
    for i in range(nbe):
        s1 = coord[int(tri[i][0])]
        s2 = coord[int(tri[i][1])]
        s3 = coord[int(tri[i][2])]
        print("s1:",s1)
        [pas,qualite] = pas_et_qualite_triangle(s1, s2, s3)
        if pas > maxPas :
            maxPas = pas
        if maxQual < qualite :
            maxQual = qualite
    print(maxPas,maxQual)
    return maxPas,maxQual
def mes_t(s1,s2,s3) :
    aire = 0.5*abs(cross(s2-s1,s3-s1))
    return aire

def mes_a(s1,s2) :
    return linalg.norm(s1-s2)

def fct_u(x,y) :
    return 1+math.sin((math.pi/2)*x)+x*(x-4)*math.cos((math.pi/2)*y)

def fct_uE(x,y) :
    return Ue

def fct_f(x,y) :
    return (pow(math.pi,2))/4*math.sin((math.pi/2)*x)+(((pow(math.pi,2)/4)*pow(x,2))-pow(math.pi,2)*x-2)*math.cos((math.pi/2)*y)
def fct_kappa() :
    return Keau

def fct_alpha() :
    return alpha

def coeffelem_P1_rigid(s1,s2,s3) :
    retour = zeros((3,3))
    aire = mes_t(s1,s2,s3)
    k = fct_kappa()
    retour[0][0] = (k/(4*aire))*(pow(s2[0]-s3[0],2)+pow(s2[1]-s3[1],2))
    retour[0][1] = (k/(4*aire))*(-(s1[0]-s3[0])*(s2[0]-s3[0])-(s1[1]-s3[1])*(s2[1]-s3[1]))
    retour[0][2] = (k/(4*aire))*(-(s3[0]-s2[0])*(s1[0]-s2[0])-(s3[1]-s2[1])*(s1[1]-s2[1]))
    retour[1][0] = retour[0][1]
    retour[1][1] = (k/(4*aire))*(pow(s3[0]-s1[0],2)+pow(s3[1]-s1[1],2))
    retour[1][2] = (k/(4*aire))*(-(s2[0]-s1[0])*(s3[0]-s1[0])-(s2[1]-s1[1])*(s3[1]-s1[1]))
    retour[2][0] = retour[0][2]
    retour[2][1] = retour[1][2]
    retour[2][2] = (k/(4*aire))*(pow(s1[0]-s2[0],2)+pow(s1[1]-s2[1],2))
    return retour
def fonction_masse(s1,s2,s3):
    retour = zeros((3,3))
    aire = mes_t(s1,s2,s3)
    #print("aire:",aire)
    retour[0][0] = aire/3
    retour[1][1] = aire/3
    retour[2][2] = aire/3
    return retour

def coeffelem_P1_source(T) :
    vec_source = np.ones(3)
    x = zeros(3,float)
    y = zeros(3,float)
    x[0] = T[0][0]
    x[1] = T[1][0]
    x[2] = T[2][0]
    y[0] = T[0][1]
    y[1] = T[1][1]
    y[2] = T[2][1]
    vec_source = (mes_t(T[0], T[1], T[2])/3)*fct_f(sum(x)/3, sum(y)/3)*vec_source
    return vec_source


def coeffelem_P1_transf(A) :
    vec_trans = np.ones(2)
    x = zeros(2,float)
    y = zeros(2,float)
    x[0] = A[0][0]
    x[1] = A[1][0]
    y[0] = A[0][1]
    y[1] = A[1][1]
    vec_trans = (mes_a(A[0], A[1])/2)*fct_alpha()*fct_uE(sum(x)/2,sum(y)/2)*vec_trans
    return vec_trans

def coeffelem_P1_poids(A) : 
    mat_poids = np.array([[2.,1.],[1.,2.]])
    mat_poids = (mes_a(A[0], A[1])/6)*fct_alpha()*mat_poids
    return mat_poids

def assemblage_EF_P1(nbe,nba,tri,ar,refa,nbn,coord) :
    A = np.zeros((nbn,nbn))
    M = np.zeros((nbn,nbn))
    F = np.zeros(nbn,float)
    for l in range(nbe) :
        I1 = int(tri[l][0])
        I2 = int(tri[l][1])
        I3 = int(tri[l][2])
        K = coeffelem_P1_rigid(coord[I1], coord[I2], coord[I3])
        Fl = coeffelem_P1_source(np.array([coord[I1], coord[I2], coord[I3]]))
        M1 = fonction_masse(coord[I1],coord[I2],coord[I3])
        A[I1][I1] += K[0][0]
        A[I1][I2] += K[0][1]
        A[I1][I3] += K[0][2]
        F[I1] += Fl[0]
        
        A[I2][I1] += K[1][0]
        A[I2][I2] += K[1][1]
        A[I2][I3] += K[1][2]
        F[I2] += Fl[1]
        
        A[I3][I1] += K[2][0]
        A[I3][I2] += K[2][1]
        A[I3][I3] += K[2][2]
        F[I3] += Fl[2]

        M[I1][I1] += M1[0][0]
        M[I1][I2] += M1[0][1]
        M[I1][I3] += M1[0][2]
        
        M[I2][I1] += M1[1][0]
        M[I2][I2] += M1[1][1]
        M[I2][I3] += M1[1][2]
        
        M[I3][I1] += M1[2][0]
        M[I3][I2] += M1[2][1]
        M[I3][I3] += M1[2][2]
    K = A.copy()
    for a in range (nba) :
        I1 = int(ar[a][0])
        I2 = int(ar[a][1])
        
        Pa = coeffelem_P1_poids(np.array([coord[I1], coord[I2]]))
        Ea = coeffelem_P1_transf(np.array([coord[I1], coord[I2]]))
        
        A[I1][I1] += Pa[0][0]
        A[I1][I2] += Pa[0][1]
        F[I1] += Ea[0]
        A[I2][I1] += Pa[1][0]
        A[I2][I2] += Pa[1][1]
        F[I2] += Ea[1]
        
    return(A,F,K,M)

def graphes(nbn,U,Uh,X,Y) :
    #figure(2)
    axes = plt.axes(projection="3d")
    axes.set_zlabel("erreur Eh")
    axes.plot_trisurf(X,Y,U-Uh,linewidth = 0.2,cmap=plt.cm.CMRmap)
    plt.xlabel("x")
    plt.ylabel("y")
    
    plt.show()

#Algo projet : ite = nb itération , nbn = nb sommet
def Algo_iteratif(ite , nbn : int , nbe :int , nba :int ,tri, ar , coord, refa):
    vec = np.ones(nbn)
    Gold = np.zeros(nbn)
    Uold = Uc * vec
    Utemps = np.zeros((ite+1,2))
    t = 0
    #init data
    #Recuperation des paramètres A, F , K
    [A,F,K,M] = assemblage_EF_P1(nbe,nba,tri,ar,refa,nbn,coord)
    B = (c/delta)*M + A
    Utemps[0]=[0,mean(Uold)-Kelvin]

    plt.ion()

    for i in range(1,ite+1):
        t = t + delta
        d = t * V
        Gold = F + (c/delta) * M@Uold 
        Unew =(np.linalg.inv(B))@Gold
        #print("Vold:",i," = ", Vnew)
        Uold = Unew
        Utemps[i]=[d,mean(Uold)-Kelvin]

        #affichage a chaque D de la canalisation
        fig = plt.figure(1)
        plt.title(f"Température moyenne {round(mean(Uold-Kelvin),3)}")
        if (d - int(d))>0.5:
            dist = int(d)+1
        else:
            dist = int(d)
        plt.suptitle(f"Section du conduit de climatisation à {int(dist)} m")
        trisurf = plt.tripcolor(coord[:, 0], coord[:, 1], Uold-Kelvin,
        cmap=plt.get_cmap("hot_r"),
        vmin=Uc-Kelvin,
        vmax=Ue-Kelvin,
        linewidth=0.2,
        antialiased=True)
        plt.pause(1)
        #sauvegarde dans le fichier save.
        plt.savefig("./save/section"+str(i) + ".png")
        plt.clf()

        
    plt.close()
    plt.ioff()
    #affichage console du min / Moyenne / max de la temperature
    print("Max :",max(Uold)-Kelvin)
    print("Moyenne :", mean(Uold)-Kelvin)
    print("min :", min(Uold)-Kelvin)

    # Création de la colormap
    my_cmap = plt.get_cmap('hot_r')

    #affichage figure final avec graphe
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    # Creating plot
    fig.set_size_inches(15, 5)
    trisurf = ax1.tripcolor(coord[:, 0], coord[:, 1], Uold-Kelvin,
        cmap=my_cmap,
        vmin=Uc-Kelvin,
        vmax=Ue-Kelvin,
        linewidth=0.1,
        antialiased=True)
    ax1.set_title(f"Coupe de canalisation à la fin du segments \n Moyenne : {mean(Uold)-Kelvin}")
    
    fig.colorbar(trisurf, shrink=0.5, aspect=5)
    
    ax2.plot(Utemps[:, 0],Utemps[:, 1] )
    ax2.set_title("Graphe de la température moyenne dans la canalisation")
    ax2.set_xlabel("Distance parcouru en mètre")
    ax2.set_ylabel("Température moyenne")
    plt.show()
    fig.savefig("./save/sectionfinal.png")

#main lancement
if __name__ == '__main__':
   
    nbn,nbe,nba,coord,tri,ar,refn,reft,refa = lit_fichier_msh(FichierMaillage)
    Algo_iteratif(n,nbn,nbe,nba,tri,ar,coord,refa)