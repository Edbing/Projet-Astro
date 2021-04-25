import numpy as np
from matplotlib import pyplot as plt

mydir = '/Users/edgarbingler/Desktop/Projet tutoré/Raw data/Isochrones NGC 6752/'
myfile = 'PadovaIsochrones_for_ngc6752_HST_vegamags.dat'

isoch_all = np.loadtxt(mydir+myfile, unpack=True)

with open(mydir+'Isoch_titles.txt') as ff :
    isoch_titles = np.array(  ff.read().splitlines() )
    # I have converted to an np.array because I like using np.where

logTe='logTe'
logg='logg'
F438W='F438Wmag'
F336W='F336Wmag'
F275W='F275W1mag'
F814W='F814Wmag'
idx1 = np.where(isoch_titles==logTe)[0][0]
idx2 = np.where(isoch_titles==logg)[0][0]
idx438 = np.where(isoch_titles==F438W)[0][0]
idx336 = np.where(isoch_titles==F336W)[0][0]
idx275 = np.where(isoch_titles==F275W)[0][0]
idx814 = np.where(isoch_titles==F814W)[0][0]
idxmet = np.where(isoch_titles=='MH')[0][0]
idxage = np.where(isoch_titles=='logAge')[0][0]
idxmini = np.where(isoch_titles=='Mini')[0][0]
jdxsolar = np.where( (isoch_all[idxmet,:]>-0.1) & (isoch_all[idxage,:]==10.) )[0]
#jdxsubsol = np.where(isoch_all[idxmet,:]<-0.1) 

distancemodulus = 12  # = 5 * log10 (distance in units of 10 pc)

plt.figure(2); plt.clf()
#plt.scatter(isoch_all[idx2,:],isoch_all[idx2,:])
plt.scatter(np.exp((isoch_all[idx1,jdxsolar])*np.log(10)),np.exp((isoch_all[idx2,jdxsolar])*np.log(10)), s=0.8)
#plt.scatter(isoch_all[idx1,jdxsubsol], isoch_all[idx1,jdxsubsol]+distancemodulus, s=0.3)
plt.grid()
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.xlabel('Teff')
plt.ylabel('g')
plt.show()

esp = 100
plt.figure(2); plt.clf()
#plt.scatter(isoch_all[idx2,:],isoch_all[idx2,:])
plt.scatter(np.exp((isoch_all[idx1,jdxsolar])*np.log(10))/esp,isoch_all[idx2,jdxsolar], s=0.8,)
#plt.scatter(isoch_all[idx1,jdxsubsol], isoch_all[idx1,jdxsubsol]+distancemodulus, s=0.3)
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.grid()
plt.gca().set_xticks(np.linspace(2000/esp,6000/esp,41))
plt.gca().set_yticks(np.linspace(-1,6,15))
plt.gca().set_xlim(6000/esp,2000/esp)
plt.xlabel('Teff *10^-2')
plt.ylabel('log g')
plt.show()

plt.figure(2); plt.clf()
#plt.scatter(isoch_all[idx2,:],isoch_all[idx2,:])
plt.scatter((isoch_all[idxmini,jdxsolar]),isoch_all[idx2,jdxsolar], s=0.8,)
#plt.scatter(isoch_all[idx1,jdxsubsol], isoch_all[idx1,jdxsubsol]+distancemodulus, s=0.3)
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.grid()
plt.gca().set_xticks(np.linspace(0.07,1.05,15))
plt.gca().set_yticks(np.linspace(-1,6,15))
#plt.gca().set_xlim(6000/esp,2000/esp)
plt.xlabel('Masse initiale')
plt.ylabel('log g')
plt.show()

plt.scatter((isoch_all[idx275,jdxsolar]-isoch_all[idx336,jdxsolar])-(isoch_all[idx336,jdxsolar]-isoch_all[idx438,jdxsolar]),isoch_all[idx336,jdxsolar])
plt.grid()

plt.ylabel('Mag F336W')
plt.xlabel('Color (F275W-F336W)-(F336W-F438)')
plt.show()

plt.scatter(isoch_all[idx438,jdxsolar]-isoch_all[idx814,jdxsolar],isoch_all[idx814,jdxsolar])
plt.plot(isoch_all[idx438,jdxsolar]-isoch_all[idx814,jdxsolar],isoch_all[idx814,jdxsolar])
plt.grid()

plt.ylabel('Mag F814W')
plt.xlabel('Color (F438W-F814W)')
plt.gca()
plt.show()
print(isoch_titles)