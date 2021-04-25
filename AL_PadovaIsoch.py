import numpy as np
from matplotlib import pyplot as plt

mydir = '/Users/edgarbingler/Desktop/Projet tutoré/Raw data/Isochrones NGC 6752/'
myfile = 'PadovaIsochrones_for_ngc6752_HST_vegamags.dat'

isoch_all = np.loadtxt(mydir+myfile, unpack=True)

with open(mydir+'Isoch_titles.txt') as ff :
    isoch_titles = np.array(  ff.read().splitlines() )
    # I have converted to an np.array because I like using np.where

magred='F438Wmag'
magblue='F275W1mag'
idx1 = np.where(isoch_titles==magred)
idx2 = np.where(isoch_titles==magblue)
idxmet = np.where(isoch_titles=='MH')
jdxsolar = np.where(isoch_all[idxmet,:]>-0.1)
jdxsubsol = np.where(isoch_all[idxmet,:]<-0.1) 

distancemodulus = 12  # = 5 * log10 (distance in units of 10 pc)

plt.figure(2); plt.clf()
#plt.scatter(isoch_all[idx2,:],isoch_all[idx2,:])
plt.scatter(isoch_all[idx2,jdxsolar]-isoch_all[idx1,jdxsolar], isoch_all[idx1,jdxsolar]+distancemodulus, s=0.4)
plt.scatter(isoch_all[idx2,jdxsubsol]-isoch_all[idx1,jdxsubsol], isoch_all[idx1,jdxsubsol]+distancemodulus, s=0.3)
plt.grid()
plt.gca().invert_yaxis()
plt.xlabel('color')
plt.ylabel('magnitude')
plt.show()

print(isoch_titles)