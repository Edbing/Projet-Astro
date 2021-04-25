#! /Users/lancon/anaconda3/bin/python
#
# Manipulation des fichiers hlsp_hugs_hst_...._ngc6752_....txt
# (qui sont des fichiers de type texte, c'est-a-dire ascii).


# Imports
import numpy as np
from matplotlib import pyplot as plt

# Basic functions
def read_catalog(myfile,mytitles) :
   data = np.loadtxt(myfile, unpack=False, comments=['#','R'])
   # NB: the 'R' in the above means that in every line you stop reading
   #     when you encounter character 'R'. This is because the star 
   #     identifications in the second-to-last column start with R and 
   #     loadtxt wants only numbers. We could use another reading method,
   #     but we don't need the star IDs anyway, we can do without the 
   #     last column as well. 
   #     Because we reject the last two columns, I do the same for the
   #     column titles.
   with open(mytitles) as ff :
       titles = ( ff.read().splitlines() )[0:-2]
   return titles, data

# Directories and files
mydir='/Users/edgarbingler/Desktop/Projet tutoré/Raw data/NGC 6752/'
file2='hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc6752_multi_v1_catalog-meth2.txt'
mytitles='colnames_meth2.txt'
#  colnames_meth2.txt was made from hlsp...-meth2.txt by hand, in a text
#  editor, and contains just one column name per line, in the order that the
#  columns appear in file2.

# Main program :
#---------------

titles, data = read_catalog(mydir+file2, mydir+mytitles)
#     The above gives you data, with shape(data)=(54555, 35),
#     and titles, with shape(titles)=(35,)

# Associate the column names with a column number
idxc = { i : j for j,i in enumerate(titles) }
print(idxc)
# Now idxc['F275W_mag'] is the column number for magnitude F275W
mask1 = np.logical_and(np.array(data[:,idxc['F438W_mag']]) >= 5, np.array(data[:,idxc['F438W_mag']]) <= 35)
mask2 = np.logical_and(np.array(data[:,idxc['F438W_mag']]) >= 5, np.array(data[:,idxc['F438W_mag']]) <= 35)
# If you wish to work with numpy arrays, make the ones you need :
F275W_mag = np.array(data[:,idxc['F275W_mag']])
F814W_mag = np.array(data[:,idxc['F814W_mag']])
F275W_err = np.array(data[:,idxc['F275W_RMS']]) # 1-sigma uncertainty on that mag
F438W_mag = np.array(data[:,idxc['F438W_mag']])
F275WmF438W = F438W_mag - F814W_mag
#.... make all the arrows you need
# (instead of this, you could just set npdata = np.array(data) 
# and access columns of that array directly for plots or calculations,
# without extracting 1-dimensional arrays first... )
mask1 = np.logical_and(np.array(F438W_mag) >= 5, F438W_mag <= 35)
mask2 = np.logical_and(F275WmF438W >= -10, F275WmF438W <= 10)
# Make some plots 
plt.figure(0); plt.clf()

np.shape(F275W_mag)

ax = plt.scatter(F438W_mag-F814W_mag, F438W_mag, alpha=0.3, s=0.3)
plt.xlabel('F438W-F814W')
plt.ylabel('F438W')
plt.gca().set_ylim(10,32.5)
plt.gca().set_xlim(-1,6)
plt.gca().invert_yaxis()

plt.show()
# You will notice lots of points have absurd values. The interesting 
# values for the colors are typically between -10 and 10, 
# and for the mags between 5 and 35.
# At high magnitudes (low brightness) there is a "cloud" of points that 
# are not interesting because they are too "noisy" (small signal-to-noise ratio)
# You can make a mask to select the values you want to use, or you
# can just set the y-axis and x-axis of the plots so the outliers are
# outside the figure (without being removed from your data arrays).




