 #This file try to analysis those data from Scatter_data_1 in D:/mask/MUSIC_redshift
import sys
sys.path.insert(0,'D:/mask')
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import glob
import h5py 
import matplotlib.pyplot as plt
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
#read the value of redshift
f = h5py.File('Redshift.h5','r')
y0 = f['a']
z = np.array(y0)
f.close()
print(z)
#read the main_tree text(which saved in the Merg_tree_read file)
with h5py.File('main_tree.h5') as f:
    y = f['a']
    main_tree = np.array(y)
print(main_tree[0,:])