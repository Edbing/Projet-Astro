
import astropy.io.fits as pyfits #Permet de manipluer des .fits
import numpy as np  
import matplotlib.pyplot as plt # Pour tracer les graphes      
from scipy import signal as sig # Pour lisser les graphes      
import math # Pour le calcul
from scipy.optimize import leastsq

VISdir = r'/Users/edgarbingler/Desktop/Projet tutoré/Raw data/Z-0.0/' # Raccourci l'écriture
Teff = '05900' # Température effective voulue en Kelvin
grav1 = '-0.50' # Gravité de surface voulue en cm.s^-2 (à noter que le signe est inversé)
grav2 = '-2.00'
comp = '-0.0' # Composition de l'étoile ([Fe/H] ou [M/H])
comp2 = '-1.5'

def lorentzian( x, x0, a, gam ):
    return a * gam**2 / ( gam**2 + ( x - x0 )**2)

def fitsfile(x,y,z):
    return VISdir+'lte'+x+y+z+'.fits' # On crée une fonction car c'est plus pratique

def hh(x,y,z):
    return pyfits.open(fitsfile(x,y,z)) #On ouvre le fichier

def hdu(x,y,z):
    return hh(x,y,z)[0]

def data(x,y,z):
    return hdu(x,y,z).data

#hh(Teff,grav1,comp).info() # Donne des infos sur le fichier .fits

toto =hdu(Teff,grav1,comp).data.shape

print(hdu(Teff,grav1,comp).header)

def fluxbrut(x,y,z):
    return data(x,y,z)/np.median(data(x,y,z)) # Dans le .fits on a les infos de flux en fonction de lambda
                                # On normalise pour rendre le graphe plus agréable
def thistype(x,y,z):
    return hdu(x,y,z).header.get('CTYPE1')                               

#print("Wavelength scale in fits header : ", thistype(Teff,grav1,comp))

def lambda0(x,y,z):
    return hdu(x,y,z).header.get('CRVAL1')

def dlambda(x,y,z): 
    return hdu(x,y,z).header.get('CDELT1')

def npix(x,y,z):
    return hdu(x,y,z).header.get('NAXIS1')

def lambdaend_plus1(x,y,z):
    return lambda0(x,y,z) +npix(x,y,z)*dlambda(x,y,z) # On crée la variable qui correspond au dernier lambda voulu

def waves(x,y,z):
    return np.exp(np.arange(lambda0(x,y,z),lambdaend_plus1(x,y,z),dlambda(x,y,z)))

#print(waves(Teff,grav1,comp)[0:5], waves(Teff,grav1,comp)[-5:])   # 5 first ones and 5 last ones

def fluxlisse(x,y,z):
    
    kernel = sig.gaussian(100,10.)/math.sqrt(2.*math.pi*100.) # On définie la fonction pour le lissage
    #kernel = lorentzian( 0, 0,10, 10 )/math.sqrt(2.*math.pi*100.) # On définie la fonction pour le lissage
    return sig.fftconvolve(fluxbrut(x,y,z),kernel,mode='same') # On lisse le flux

def tracageg(x,y,z,col,off,dec):
    plt.plot(waves(x,y,z)*1.000200, fluxbrut(x,y,z)+off, 'k', linewidth=1.0)      
    plt.plot(waves(x,y,z)*1.000200, fluxlisse(x,y,z)+off, col, linewidth=1.5, label='g = '+y)
    plt.ylim(-0.1,1.5*np.percentile(fluxbrut(x,y,z),99))  # On ne prends pas en compte le dernier pourcentile pour exclure de potentielles abérations
    plt.xlim(2000,25000-dec)
    plt.xlabel("longueur d'onde",fontsize=25)
    plt.ylabel("Intensité lumineuse",fontsize=25)

def tracaget(x,y,z,col,off):
    plt.plot((waves(x,y,z)*1.000200), fluxbrut(x,y,z)+off, 'k', linewidth=1.0)      
    plt.plot((waves(x,y,z)*1.000200), fluxlisse(x,y,z)+off, col, linewidth=1.5, label='Teff = '+x)
    plt.ylim(-0.1,1.5*np.percentile(fluxbrut(x,y,z),99))  # On ne prends pas en compte le dernier pourcentile pour exclure de potentiels abérations
    plt.xlabel("longueur d'onde",fontsize=25)
    plt.ylabel("Intensité lumineuse",fontsize=25)
    
def tracagec(x,y,z,col,off):
    plt.plot(waves(x,y,z)*1.000200, fluxbrut(x,y,z)+off, 'k', linewidth=1.0)      
    plt.plot(waves(x,y,z)*1.000200, fluxlisse(x,y,z)+off, col, linewidth=1.5, label='[M/H] = '+z)
    plt.ylim(-0.1,1.5*np.percentile(fluxbrut(x,y,z),99))  # On ne prends pas en compte le dernier pourcentile pour exclure de potentiels abérations
    plt.xlabel("longueur d'onde",fontsize=25)
    plt.ylabel("Intensité lumineuse",fontsize=25)
        
plt.clf()
plt.figure(figsize=(15,15))             
tracageg(Teff,grav1,comp,'r',0,0)
tracageg(Teff,'-3.00',comp,'g',5,0)
tracageg(Teff,'-6.00',comp,'b',10,0)
plt.title("Différence de spectre en fonction de la gravité de surface à une température (en K) de "+Teff,fontsize=28)
plt.gca().set_ylim(0, 14)
plt.legend(prop={'size':25})
plt.grid()         
plt.show()
#plt.savefig(fname='Infgrav',format='png')

plt.figure(figsize=(15,15))             
tracageg('10400','-2.00',comp,'r',0,8000)
tracageg('10400','-4.00',comp,'g',9,8000)
tracageg('10400','-6.00',comp,'b',18,8000)
plt.title("Différence de spectre en fonction de la gravité de surface à une température (en K) de "+Teff,fontsize=28)
plt.gca().set_ylim(0, 28)
plt.legend(prop={'size':25})
plt.grid()         
plt.show()

plt.figure(figsize=(15,15)) 
tracaget(Teff,grav2,comp,'r',0)
tracaget('07200',grav2,comp,'g',4)
tracaget('10400',grav2,comp,'b',9)
plt.title("Différence de spectre en fonction de la température effective",fontsize=28)
plt.gca().set_ylim(0, 20)
plt.legend(prop={'size':25})
plt.grid()         
plt.show()
#plt.savefig(fname='InfTeff.pdf')

plt.figure(figsize=(15,15)) 
tracagec(Teff,grav1,comp,'r',0)
VISdir = r'/Users/edgarbingler/Desktop/Projet tutoré/Raw data/Z-1.5/'
tracagec(Teff,grav1,comp2,'g',5)
plt.title("Différence de spectre en fonction de la composition à une température (en K) de "+Teff,fontsize=28)
plt.gca().set_ylim(0, 9)
plt.legend(prop={'size':25})
plt.grid()         
plt.show()
#plt.savefig(fname='Infcomp',format='png')

VISdir = r'/Users/edgarbingler/Desktop/Projet tutoré/Raw data/Z-0.0/'
plt.figure(figsize=(15,15)) 
tracagec('10400',grav2,comp,'r',0)
VISdir = r'/Users/edgarbingler/Desktop/Projet tutoré/Raw data/Z-1.5/'
tracagec('10400',grav2,comp2,'g',9)
plt.title("Différence de spectre en fonction de la composition à une température de 10400K",fontsize=28)
plt.gca().set_ylim(0, 20)
plt.legend(prop={'size':25})
plt.grid()         
plt.show()

