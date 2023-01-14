import pylab
from pylab import *
from PIL import *
from PIL import Image
import numpy as np
from math import *

def sobel(imgtab : np.array) -> np.array :
    tab = np.copy(imgtab)
    for i in range(imgtab.shape[0]) :
        for j in range(imgtab.shape[1]) :
            Gx = Gy = 0
            
            # Calcul de Gx
            if(i-1 >= 0 and j-1 >= 0 and j+1 < imgtab.shape[1]) :
                Gx += -1 * int(imgtab[i-1][j-1]) + 1 * int(imgtab[i-1][j+1])
            if(j-1 >= 0 and j+1 < imgtab.shape[1]) :
                Gx += -2 * int(imgtab[i][j-1])   + 2 * int(imgtab[i][j+1])
            if(i+1 < imgtab.shape[0] and j-1 >= 0 and j+1 < imgtab.shape[1]) :
                Gx += -1 * int(imgtab[i+1][j-1]) + 1 * int(imgtab[i+1][j+1])
            
            # Calcul de Gy
            if(i-1 >= 0 and j-1 >= 0 and j+1 < imgtab.shape[1]) :
                Gy += -1 * int(imgtab[i-1][j-1]) - 2 * int(imgtab[i-1][j]) - 1 * int(imgtab[i-1][j+1])
            if(i+1 < imgtab.shape[0] and j-1 >= 0 and j+1 < imgtab.shape[1]) :
                Gy += 1  * int(imgtab[i+1][j-1]) + 2 * int(imgtab[i+1][j]) + 1 * int(imgtab[i+1][j+1])
            
            # Calcul de G
            G = sqrt(Gx**2 + Gy**2)
            
            tab[i][j] = G
    return tab

def detectionPCC(imgtabsobel : np.array,nbsupp : int , img : Image) -> Image:
    hauteur = imgtabsobel.shape[0]-1
    largeur = imgtabsobel.shape[1]
    poid= np.zeros ((largeur))
    chemin = np.zeros((largeur,hauteur+1))
    
    #i = indice du pixel de base
    for i in range(0,largeur):
        ind = i
        poid[ind] = 0
        poid[ind] += imgtabsobel[hauteur][ind]
        chemin[i][hauteur] = ind
        for h in range(hauteur-1,-1,-1):
            pix = [1000,1000,1000]
            #verification que l'on deborde pas vers le haut h et i vers la gauche
            if(ind-1>=0 and h>0):
                pix[0] = imgtabsobel[h][ind-1]
            
            if(ind>=0 and h>0):
                pix[1] = imgtabsobel[h][ind] 
            
            if(ind+1<largeur and h>0):
                pix[2] = imgtabsobel[h][ind+1] 
            
            poid[i] += min(pix)  
            if(min(pix) == pix[0]) : 
                ind = ind-1
            if(min(pix) == pix[2]) : 
                ind = ind+1
            chemin[i][h] = ind

    for suppression in range(nbsupp):
        indiceasupp = indPCC(poid)
        newpoid= np.delete(poid,indiceasupp,0)
        poid = np.copy(newpoid)
        img = delPCC(img,chemin,indiceasupp)
    return img


def delPCC(img : Image , chemin: np.array, indiceasupp : int) -> np.array:
    hauteur = img.height
    largeur = img.width
    newimg =  img.copy()
    for h in range(hauteur-1):
        #chemin[indiceasupp][h] renvoie l indice a supp sur la ligne
        #newimgtab[h] = decaleligne(img,int(chemin[indiceasupp][h]),largeur-1)
        for y in range(int(chemin[indiceasupp][h]) , largeur-1):
             if(y < largeur):
                newimg.putpixel((y,h),img.getpixel((y+1,h))) 
       
    return newimg

def resize(img : Image, nbsupp) -> Image:
    hauteur = img.height
    largeur = img.width
    img1 = img.crop((0,0,largeur-nbsupp,hauteur))

    return img1          
    
#donne l indice du plus court chemin dans le tableau / dans l image
def indPCC(poid : np.array) -> int:
    for i in range(poid.shape[0]):
        if poid[i] == min(poid):
            return i

    
if __name__ == "__main__":
    chemin = "tower_480.png"
    nbsuppcol = 75
    nbsuppligne = 20
    #Recuperation de l image
    img = Image.open(chemin)
    newimg = img.copy()
    if(nbsuppcol>0):
        tabimg = np.array(img)
        #conversion en nuance de gris
        tabimggris = np.array(img.convert('L'))
        #Image avec filtre sobel
        tabsobel = sobel(tabimggris)
        Image.fromarray(tabsobel).show()
        Image.fromarray(tabsobel).save("img_sobel.png")
        #detection plus court chemin et decalage des pixels
        newimg = detectionPCC(tabsobel,nbsuppcol,img)
        #affichage image final
        newimg = resize(newimg,nbsuppcol)
    if(nbsuppligne>0):
        newimg = newimg.rotate(90,expand = True)
        #conversion en nuance de gris
        tabimggris = np.array(newimg.convert('L'))
        #Image avec filtre sobel
        tabsobel = sobel(tabimggris)
        #detection plus court chemin et decalage des pixels
        newimg = detectionPCC(tabsobel,nbsuppligne,newimg)
        #affichage image final
        newimg = resize(newimg,nbsuppligne)

        newimg = newimg.rotate(-90,expand = True)
    newimg.show()
    newimg.save("img_seamcarving.png")
