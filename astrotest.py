#! /Users/lancon/anaconda3/bin/python

# This quick-and-dirty program shows how to produce basic plots
# to display theoretical spectra (from Husser et al. 2013, version v2)
# and empirical spectra (from XSL iDR2).
#
# HISTORY : created by A.L., 2016 then modified 2018.
#
# First improvements needed : 
#    - make a python function to read an XSL spectrum (given the filename)
#      and return flux, waves, dlambda, etc
#    - make a python function to read a model spectrum (given the filename)
#      and return fluxmod, wavemod, etc
#    - make a python function that makes a certain type of plot,
#      (using the previous functions to read data).
# Then decide what figures exactly you want to make, produce them (with 
# a program, not by hand) and save them.


#  import pyfits   
#   [version plus moderne : import astropy.io.fits as fits, ou as pyfits]
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
fitsfile = VISdir+'lte02300+0.50-0.0.fits' #lte"tempeff""gravité""composition"

hh = pyfits.open(fitsfile)
hh.info()  #   [gives you info on the organization of the contents in HDUs]

hdu=hh[0] #   [Necessary to go to the first HDU]
          #   [That first HDU contains the spectrum.]
hdu.data.shape     #   [Gives you the format of the data. 
                   #    hdu.data is the data, and behaves like a Numpy array.]
                   #   Alternative without hdu=hh[0] :   hh[0].data.shape
hdu.data[0:10]  # gives you the first 10 values of data (0 to 9)
                # These are arbitrary units of energy/time/wavelength/area
                # but typically of order 1.e-13
flux = hdu.data/np.median(hdu.data)

hdu.header.items()  # gives you the header as a list of dictionaries (?):
       #  For instance : hdu.header.items()[3] is :   ('NAXIS1', 24750)
# hdu.header.has_key('NAXIS1') #  gives you True.
hdu.header.get('NAXIS1')      # returns the value 24750
       # Other syntax : hdu.header['NAXIS1']
# hdu.header.values()          #  gives you the header values without any comments.
       # e.g.:   hdu.header.values()[3] is 24750  (= value of key NAXIS1)

#Pour reconstruire les longueurs d'ondes a partir de NAXIS1, 
#CRPIX1, CRVAL1, CDELT1 (=CD1_1), CTYPE1 (='LINEAR' en general mais pas tjrs) :
thistype = hdu.header.get('CTYPE1')
print("Wavelength scale in fits header : ", thistype)
lambda0 = hdu.header.get('CRVAL1')
dlambda = hdu.header.get('CDELT1')
npix = hdu.header.get('NAXIS1')
lambdaend_plus1= lambda0 +npix*dlambda  # _plus1 a cause du fonctionnement de arange
#print(lambda0)
#print(lambdaend_plus1)
#print(dlambda)
waves = np.exp(np.arange(lambda0,lambdaend_plus1,dlambda))    # nm
print(waves[0:5], waves[-5:])   # 5 first ones and 5 last ones

# The convolution below is done with scipy.signal.fftconvolve
kernel = sig.gaussian(100,10.)/math.sqrt(2.*math.pi*100.)
flux2 = sig.fftconvolve(flux,kernel,mode='same')  # bugs if flux had not been normalized
# Alternative (better?) : replace with the convolution options of astropy.convolution

plt.clf()             # Clears the graphics window if it is still present
plt.plot(waves*1.000200, flux, 'k', linewidth=1.0)      # Small velocity (Doppler) shift
plt.plot(waves*1.000200, flux2, 'r', linewidth=1.5)
          # or for all lines: plt.setp(lines, linewidth=1.0)
plt.ylim(-0.1,1.5*np.percentile(flux,99))  # Adjust y axis to avoid showing outliers
          # or : plt.axis([xmin,xmax,ymin,ymax])
plt.show()
#plt.save #penser à faire des grandes légendes directement pour éviter les problèmes de scaling

# O.Husser model spectra :
#==========================
#OHdir = '/Work/lancon/data/Husser/Z-0.0/'  
#fitsmod = OHdir+'lte03600-0.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
OHdir = r'C:\\Users\\edgar\\Desktop\\Projet tutoré\\Raw data\\Z-0.0\\'
fitsmod = OHdir+'lte05200-3.00-0.0.fits'
hhmod = pyfits.open(fitsmod)
hdumod = hhmod[0]   # Goes to first (and only) extension of these files.
fluxmod = hdumod.data/np.median(hdumod.data)  # 212027 points !! 
          # Orig units are erg/s/cm2/cm and typically of order 1.e13
          # In these files, wavelengths are log-spaced : 
          # hdumod.header.get('CTYPE1') is 'AWAV-LOG'
wavemod0 = hdumod.header.get('CRVAL1')
dlogwave = hdumod.header.get('CDELT1')
npixmod = fluxmod.size
wavemod = np.exp( np.arange(wavemod0, wavemod0+npixmod*dlogwave, dlogwave) )/10. # AA -> nm 

# Smooth model
kernel = sig.gaussian(100,12.)/math.sqrt(2.*math.pi*12*12)
fluxmod2 = sig.fftconvolve(fluxmod,kernel,mode='same')

# Normalize model to value of data near 0.85mu 
f1 = flux[ (waves>=800)*(waves<=900) ].mean()
f2 = fluxmod2[ (wavemod>=800)*(wavemod<=900) ].mean()


plt.plot(wavemod, fluxmod*f1/f2, '--', linewidth=1.0, color='slategrey')
plt.plot(wavemod, fluxmod2*f1/f2, '--', linewidth=1.5, color='salmon')
plt.xlim(500,1100)  
#plt.show()


# Other figure, with two axis systems that share coordinates :
#  (cf matplotlib.org/examples/pylab_examples/subplots_demo.html)
thisframe, (ax1, ax2) = plt.subplots(2, 1, sharey=True, sharex=True, figsize=(5,7))
ax1.plot(waves*1.0002, flux, 'k', linewidth=1.0)
ax1.plot(wavemod, fluxmod*f1/f2, '-', linewidth=1.0, color='slategrey')
ax1.set_ylabel('Original')
ax2.plot(waves*1.0002, flux2, 'r', linewidth=1.0)
ax2.plot(wavemod, fluxmod2*f1/f2, '-', linewidth=1.0, color='salmon')
ax2.set_ylabel('Smoothed')
ax2.set_xlabel('Wavelength')
plt.ylim(-0.1,1.5*np.percentile(flux,99))
plt.xlim(500,1100) 
plt.show()

# Save figure :  (saves the last figure plotted or activated interactively)
#===============
#plt.savefig("ALtest.png") # ou "ALtest.pdf", ou "ALtest.eps", ou "ALtest.ps"
