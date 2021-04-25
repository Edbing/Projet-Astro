import astropy.io.fits as pyfits #Permet de manipluer des .fits
import numpy as np  
import matplotlib.pyplot as plt # Pour tracer les graphes      
from scipy import signal as sig # Pour lisser les graphes    
from scipy.interpolate import interp1d   
import math # Pour le calcul
#from scipy import integrate as intg


VISdir = r'/Users/edgarbingler/Desktop/Projet tutoré/Raw data/Z-0.0/' # Raccourci l'écriture
Teff = '05700' # Température effective voulue en Kelvin
grav1 = '-0.50' # Gravité de surface voulue en cm.s^-2 (à noter que le signe est inversé)
grav2 = '-3.00'
comp = '-0.0' # Composition de l'étoile ([Fe/H] ou [M/H])
comp2 = '-1.5'
c = 299792458 #m.s-1
Jy = 1.e23
Rsol = 7*10^10  #cm


def fitsfile(x,y,z):
    return VISdir+'lte'+x+y+z+'.fits' # On crée une fonction car c'est plus pratique

def hh(x,y,z):
    return pyfits.open(fitsfile(x,y,z)) #On ouvre le fichier

def hdu(x,y,z):
    return hh(x,y,z)[0]

def rayon(x,y,z):
    #print(hdu(x,y,z).header.get('PHXREFF'))
    return hdu(x,y,z).header.get('PHXREFF')

def fluxbrut(x,y,z):
    return hdu(x,y,z).data*((rayon(x,y,z)/Rsol)**2)*4*np.pi # Dans le .fits on a les infos de flux en fonction de lambda 
                                # On normalise pour rendre le graphe plus agréable
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

def fluxlisse(x,y,z):
    
    kernel = sig.gaussian(100,10.)/math.sqrt(2.*math.pi*100.) # On définie la fonction pour le lissage
    #return sig.fftconvolve(fluxbrut(x,y,z),kernel,mode='same') # On lisse le flux
    return fluxbrut(x,y,z)
def fluxnu(x,y,z):
    return ((fluxlisse(x,y,z)*waves(x,y,z)*waves(x,y,z))/c) 

def nu(x,y,z):
    return (c/waves(x,y,z))


def filtrew(a,b,x,y,z):
    mask = np.logical_and(waves(x,y,z) >= a, waves(x,y,z) <= b)
    fluxfiltre=fluxbrut(x,y,z)[mask]
    wavefiltre=waves(x,y,z)[mask]
    fil= wavefiltre/wavefiltre
    
    return wavefiltre,fluxfiltre,fil

def filtrenu(a,b,x,y,z):
    mask = np.logical_and(waves(x,y,z) >= a, waves(x,y,z) <= b)
    fluxnufiltre=fluxnu(x,y,z)[mask]
    nufiltre=nu(x,y,z)[mask]
    fil= nufiltre/nufiltre

    return nufiltre,fluxnufiltre,fil

 
def fnumoy(a,b,x,y,z):
   x1 = np.trapz((filtrenu(a,b,x,y,z)[1])/filtrenu(a,b,x,y,z)[0],filtrenu(a,b,x,y,z)[0])
   x2 = np.trapz((filtrenu(a,b,x,y,z)[2])/filtrenu(a,b,x,y,z)[0],filtrenu(a,b,x,y,z)[0])
   return (x1/x2)

def mag(a,b,x,y,z):
    return -2.5*np.log(fnumoy(a,b,x,y,z))
#filtre1=filtre(a,b,x,y,z)

def indice(a,b,c,d,x,y,z):
    return 2.5*np.log(fnumoy(c,d,x,y,z)/fnumoy(a,b,x,y,z))

def diagbmv():
    donmag = []
    donind = []
    #donmag.append(mag(4000,4800,'03000',grav2,comp))
    #donmag.append(mag(4000,4800,'03200',grav2,comp))
    #donmag.append(mag(4000,4800,'03500',grav2,comp))
    donmag.append(mag(4000,4800,'04000',grav2,comp))
    donmag.append(mag(4000,4800,'04500',grav2,comp))
    donmag.append(mag(4000,4800,'05000',grav2,comp))
    donmag.append(mag(4000,4800,'06000',grav2,comp))
    donmag.append(mag(4000,4800,'07000',grav2,comp))
    donmag.append(mag(4000,4800,'10000',grav2,comp))
    donmag.append(mag(4000,4800,'12000',grav2,comp))
    donmag.append(mag(4000,4800,'15000',grav2,comp))
    #donind.append(indice(4000,4800,5000,6000,'03000',grav2,comp))
    #donind.append(indice(4000,4800,5000,6000,'03200',grav2,comp))
    #donind.append(indice(4000,4800,5000,6000,'03500',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'04000',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'04500',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'05000',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'06000',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'07000',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'10000',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'12000',grav2,comp))
    donind.append(indice(4000,4800,5000,6000,'15000',grav2,comp))
    #donmag.reverse()
    #donind.reverse()
    print(donmag,donind)
    return donmag,donind

def diagHST0():
    donmag = []
    donind = []
    #donmag.append(mag(4000,4800,'03000',grav2,comp))
    #donmag.append(mag(4000,4800,'03200',grav2,comp))
    #donmag.append(mag(4000,4800,'03500',grav2,comp))
    donmag.append(fil475('04000',grav2,comp))
    donmag.append(fil475('04500',grav2,comp))
    donmag.append(fil475('05000',grav2,comp))
    donmag.append(fil475('06000',grav2,comp))
    donmag.append(fil475('07000',grav2,comp))
    donmag.append(fil475('10000',grav2,comp))
    donmag.append(fil475('12000',grav2,comp))
    donmag.append(fil475('15000',grav2,comp))
    #donind.append(indice(4000,4800,5000,6000,'03000',grav2,comp))
    #donind.append(indice(4000,4800,5000,6000,'03200',grav2,comp))
    #donind.append(indice(4000,4800,5000,6000,'03500',grav2,comp))
    donind.append(indice1('04000',grav2,comp))
    donind.append(indice1('04500',grav2,comp))
    donind.append(indice1('05000',grav2,comp))
    donind.append(indice1('06000',grav2,comp))
    donind.append(indice1('07000',grav2,comp))
    donind.append(indice1('10000',grav2,comp))
    donind.append(indice1('12000',grav2,comp))
    donind.append(indice1('15000',grav2,comp))
    return donmag,donind


def diaglogLT():
    donT = []
    donlum = []
    donT.append(hdu('03000',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('03200',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('03500',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('04000',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('04500',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('05000',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('06000',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('07000',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('10000',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('12000',grav2,comp).header.get('PHXTEFF'))
    donT.append(hdu('15000',grav2,comp).header.get('PHXTEFF'))
    donlum.append(hdu('03000',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('03200',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('03500',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('04000',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('04500',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('05000',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('06000',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('07000',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('10000',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('12000',grav2,comp).header.get('PHXLUM'))
    donlum.append(hdu('15000',grav2,comp).header.get('PHXLUM'))
    
    #donmag.reverse()
    #donind.reverse()
    return donT,donlum
    

def openf(filt):
    
    filterW, filterT =np.loadtxt(filt, unpack=True)
    return filterT, filterW

def fil475(x,y,z):
   
    filter = "/Users/edgarbingler/Desktop/Projet tutoré/Filtres/HST_WFC3_UVIS2.F438W.dat"
    filtrans = openf(filter)[0]
    fillambd = openf(filter)[1]
    
    mask = np.logical_and(waves(x,y,z) >= fillambd[0], waves(x,y,z) <=fillambd[-1])
    Flambda = fluxlisse(x,y,z)[mask]
    lam = waves(x,y,z)[mask]
    filtre = interp1d(fillambd, filtrans)
    
    transmis = filtre(lam)
    
    
    fnu = (Flambda*lam*lam/c)* 1.e-20
    nu = (c/lam) *1.e10
    
    I1 = np.trapz((fnu*Jy/4197.72)*transmis/nu,nu) # système Vega
    I2 = np.trapz(transmis/nu,nu)
    Fluxnumoy=I1/I2
    magn = -2.5*np.log10(Fluxnumoy)
    
    #print(I1)
    return Fluxnumoy,magn,(c/lam),(Flambda*lam*lam/c),filtrans,fillambd

def fil555(x,y,z):
   
    filter = "/Users/edgarbingler/Desktop/Projet tutoré/Filtres/HST_WFC3_UVIS2.F814W.dat" #penser à rechanger
    filtrans = openf(filter)[0]
    fillambd = openf(filter)[1]
    
    mask = np.logical_and(waves(x,y,z) >= fillambd[0], waves(x,y,z) <=fillambd[-1])
    Flambda = fluxlisse(x,y,z)[mask]
    lam = waves(x,y,z)[mask]
    filtre = interp1d(fillambd, filtrans)
    
    transmis = filtre(lam)
    
    
    fnu = (Flambda*lam*lam/c)* 1.e-20
    nu = (c/lam) *1.e10
    
    I1 = np.trapz((fnu*Jy/2421.40)*transmis/nu,nu)
    I2 = np.trapz(transmis/nu,nu)
    Fluxnumoy=I1/I2
    magn = -2.5*np.log10(Fluxnumoy)
    
    
    
    return Fluxnumoy,magn,(c/lam),(Flambda*lam*lam/c),filtrans,fillambd

def indice1(x,y,z):
    return fil475(x, y, z)[1]-fil555(x, y, z)[1]

test3=fil475('05000',grav2,comp)[1]
test1=diagHST0()[1]
test2=diagHST0()[0]
    
def diagHST():
    donmag = []
    donind = []
    donmasse = []
    #donmag.append(fil475('03000',grav2,comp)[1])
    #donmag.append(fil475('03200',grav2,comp)[1])
    #donmag.append(fil475('03500',grav2,comp)[1])
    donmag.append(fil475('02800','-5.00',comp)[1])
    donmag.append(fil475('02900','-5.00',comp)[1])
    donmag.append(fil475('03000','-5.00',comp)[1])
    donmag.append(fil475('03300','-5.00',comp)[1])
    donmag.append(fil475('03700','-4.50',comp)[1])
    donmag.append(fil475('04300','-4.50',comp)[1])
    donmag.append(fil475('04700','-4.50',comp)[1])
    donmag.append(fil475('05000','-4.50',comp)[1])
    donmag.append(fil475('05200','-4.50',comp)[1])
    donmag.append(fil475('05300','-4.50',comp)[1])
    donmag.append(fil475('05400','-4.50',comp)[1])
    donmag.append(fil475('05500','-4.00',comp)[1])
    donmag.append(fil475('05600','-4.00',comp)[1])
    donmag.append(fil475('05700','-4.00',comp)[1])
    donmag.append(fil475('04900','-3.50',comp)[1])
    donmag.append(fil475('04700','-3.00',comp)[1])
    donmag.append(fil475('04600','-2.50',comp)[1])
    donmag.append(fil475('04500','-2.50',comp)[1])
    donmag.append(fil475('04300','-2.00',comp)[1])
    donmag.append(fil475('04000','-1.50',comp)[1])
    donmag.append(fil475('03800','-1.00',comp)[1])
    donmag.append(fil475('03700','-1.00',comp)[1])
    donmag.append(fil475('03500','-0.50',comp)[1])
    donmag.append(fil475('03400','-0.50',comp)[1])
    donmag.append(fil475('03200','-0.00',comp)[1])
    donmag.append(fil475('03100','-0.00',comp)[1])
    donmag.append(fil475('02800','+0.50',comp)[1])
    donmag.append(fil475('04100','-1.50',comp)[1])
    donmag.append(fil475('04400','-2.00',comp)[1])
    donmag.append(fil475('04700','-2.50',comp)[1])
    #donind.append(indice1('03000',grav2,comp))
    #donind.append(indice1('03200',grav2,comp))
    #donind.append(indice1('03500',grav2,comp))
    donind.append(indice1('02800','-5.00',comp))
    donind.append(indice1('02900','-5.00',comp))
    donind.append(indice1('03000','-5.00',comp))
    donind.append(indice1('03300','-5.00',comp))
    donind.append(indice1('03700','-4.50',comp))
    donind.append(indice1('04300','-4.50',comp))
    donind.append(indice1('04700','-4.50',comp))
    donind.append(indice1('05000','-4.50',comp))
    donind.append(indice1('05200','-4.50',comp))
    donind.append(indice1('05300','-4.50',comp))
    donind.append(indice1('05400','-4.50',comp))
    donind.append(indice1('05500','-4.00',comp))
    donind.append(indice1('05600','-4.00',comp))
    donind.append(indice1('05700','-4.00',comp))
    donind.append(indice1('04900','-3.50',comp))
    donind.append(indice1('04700','-3.00',comp))
    donind.append(indice1('04600','-2.50',comp))
    donind.append(indice1('04500','-2.50',comp))
    donind.append(indice1('04300','-2.00',comp))
    donind.append(indice1('04000','-1.50',comp))
    donind.append(indice1('03800','-1.00',comp))
    donind.append(indice1('03700','-1.00',comp))
    donind.append(indice1('03500','-0.50',comp))
    donind.append(indice1('03400','-0.50',comp))
    donind.append(indice1('03200','-0.00',comp))
    donind.append(indice1('03100','-0.00',comp))
    donind.append(indice1('02800','+0.50',comp))
    donind.append(indice1('04100','-1.50',comp))
    donind.append(indice1('04400','-2.00',comp))
    donind.append(indice1('04700','-2.50',comp))
    #donmasse.append(hdu('02800','-5.00',comp).header.get('PHXTEFF'))
    #donmag.reverse()
    #donind.reverse()
    print(donind)
    return donmag,donind

templis=['2800','2900','3000','3300','3700','4300','4700','5000','5200','5300','5400','5500','5600','5700','4900','4700','4600','4500','4300','4000','3800','3700','3500','3400','3200','3100','2800','4100','4400','4700']
gravlis=['-5','-5','-5','-5','-4.5','-4.5','-4.5','-4.5','-4.5','-4.5','-4.5','-4.5','-4','-4','-3.5','-3','-2.5','-2.5','-2','-1.5','-1','-1','-0.5','-0.5','-0','-0','+0.5','-1.5','-2','-2.5']
masslis=[0.15,0.155,0.205,0.37,0.455,0.6,0.66,0.755,0.795,0.8,0.85,0.91,0.98,1.02,1.025,1.035,1.04,1.041,1.042,1.044,1.0445,1.045,1.0455,1.0462,1.047,1.048,1.0483,1.0486,1.0489,1.049]


print(len(masslis))
"""#print(fil555('04600','-2.00','-0.0'))       
plt.figure(figsize=(15,15))
plt.plot(filtrenu(5000,6000,Teff,grav1,comp)[0],filtrenu(5000,6000,Teff,grav1,comp)[1])
#plt.plot(nu('02300',grav1,comp),fluxnu('02300',grav1,comp))
#plt.plot(nu('04600',grav2,comp),fluxnu('04600',grav2,comp))
#plt.plot(nu('07600',grav2,comp),fluxnu('07600',grav2,comp))
plt.show()
plt.plot(waves(Teff,grav1,comp),fluxbrut(Teff,grav1,comp))
plt.show()
#print(len(filtrew(5000,6000,Teff,grav1,comp)[0]))
#print(trapezes(waves,fluxbrut,2000,Teff,grav1,comp))
#integrate(5000,6000,Teff,grav1,comp)
#print(np.trapz(filtrew(12000,13000,Teff,grav1,comp)[1],filtrew(12000,13000,Teff,grav1,comp)[0]))
#print(filtrenu(5000,6000,Teff,grav1,comp)[2])
#print(mag(5000,6000,Teff,grav1,comp))"""

"""plt.plot(fil475(Teff, grav1, comp)[2],fil475(Teff, grav1, comp)[3])
plt.plot(fil475(Teff, grav1, comp)[5],fil475(Teff, grav1, comp)[4])
plt.grid()
plt.show()
plt.plot(fil555(Teff, grav1, comp)[2],fil555(Teff, grav1, comp)[3])
plt.plot(fil555(Teff, grav1, comp)[5],fil555(Teff, grav1, comp)[4])
plt.grid()
plt.show()"""

"""plt.plot(fil555(Teff, grav1, comp)[5],fil555(Teff, grav1, comp)[4],label='Filtre F555W')
plt.plot(fil475(Teff, grav1, comp)[5],fil475(Teff, grav1, comp)[4],label='Filtre F438W')
plt.xlabel("Longueur d'onde")
plt.ylabel("Coefficient de transmission")
plt.legend()
plt.grid()
plt.show()"""


#plt.plot(waves(Teff, grav1, comp),fluxlisse(Teff, grav1, comp))

"""plt.scatter(diagbmv()[1],diagbmv()[0])
plt.gca().invert_yaxis()
plt.grid()
plt.title("Diagramme couleur-magnitude d'étoiles arbitraires")
plt.ylabel("Magnitude dans le filtre F438W (unité arbitraire)",fontsize=10)
plt.xlabel("Indice de couleur mF438W - mF555W",fontsize=10)
#plt.plot(diagHST()[1],diagHST()[0])
plt.show()"""

"""plt.scatter(diagHST()[1],diagHST()[0])
plt.gca().invert_yaxis()
plt.grid()
plt.title("Diagramme Hertzsprung-Russell avec les filtres HSTF438W et HSTF814W")
plt.ylabel("Magnitude F438W (unité arbitraire)",fontsize=10)
plt.xlabel("Indice de couleur mF438W - mF814W",fontsize=10)
plt.plot(diagHST()[1],diagHST()[0])
plt.show()"""

"""plt.scatter(diagbmv()[1],diagbmv()[0],color = 'r')
plt.gca().invert_yaxis()
plt.grid()
plt.title("Diagramme Hertzsprung-Russel avec les filtres B et V")
plt.ylabel("Magnitude (unité arbitraire)",fontsize=10)
plt.xlabel("Indice de couleur B-V",fontsize=10)
plt.show()"""

"""print(hdu(Teff,grav1,comp).header)
print(fluxbrut(Teff,grav1,comp))"""


"""
plt.scatter(diagbmv()[1],diaglogLT()[0])
plt.gca().set_yscale('log')
plt.show()
plt.scatter(diagbmv()[1],diagbmv()[0])
plt.gca().invert_yaxis()
plt.grid()
plt.show()
plt.scatter(diaglogLT()[0],diaglogLT()[1])
#plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.gca().set_yscale('log')
plt.gca().set_xscale('log')
#plt.grid()
plt.show()
#plt.plot(fil475()[0],fil475()[1])

print(hdu(Teff,grav1,comp).header)"""
