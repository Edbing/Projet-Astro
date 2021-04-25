
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt       # pour les figures
from scipy import signal as sig       # pour la convolution
import math                           # pour pi, sqrt...

# XSL spectra from DR1 :
#=======================

#NIRdir = '/Volumes/ALpass1/backup_tortue_home/XSL/DR1/version_p3/NIR/'
#NIRdir = '/Users/lancon/XSL/DR1/version_p3/NIR/'
VISdir = r'C:\\Users\\edgar\\Desktop\\Projet tutoré\\Raw data\\Z-0.0\\'
#fitsfile = NIRdir+'V_AL_Mon_N_S_804309_56277.34708402_F.fits'
fitsfile = VISdir+'lte02300+0.50-0.0.fits'

hh = pyfits.open(fitsfile)
hdu =hh[0]
data = hdu.data
hh.info()
print(type(data))
print(data.shape)
hdu.data.shape
print(hdu.header[2])

#for i in range(100):
    #print(data[i])
    #i += 1""

#plt.hist(data.flat, bins=50)
#print(hdu.data[0:10])
print(hdu.header.get('NAXIS1'))
flux = hdu.data/np.median(data)

thistype = hdu.header.get('CTYPE1')
print("Wavelength scale in fits header : ", thistype)
lambda0 = hdu.header.get('CRVAL1')
dlambda = hdu.header.get('CDELT1')
npix = hdu.header.get('NAXIS1')
lambdaend_plus1= lambda0 +npix*dlambda
waves = np.arange(lambda0,lambdaend_plus1,dlambda)    # nm
print(waves[0:5], waves[-5:])   # 5 first ones and 5 last ones
kernel = sig.gaussian(100,10.)/math.sqrt(2.*math.pi*100.)
flux2 = sig.fftconvolve(flux,kernel,mode='same')

plt.clf()             
plt.plot(waves*1.000200, flux, 'k', linewidth=1.0)      
plt.plot(waves*1.000200, flux2, 'r', linewidth=1.5)
         
plt.ylim(-0.1,1.5*np.percentile(flux,99))  
plt.show()