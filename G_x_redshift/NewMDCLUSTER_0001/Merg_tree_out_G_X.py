 #This file try to analysis those data from Scatter_data_1 in D:/mask/MUSIC_redshift
import sys
sys.path.insert(0,'D:/mask')
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import glob
import h5py 
import matplotlib.pyplot as plt
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
#read the value of redshift
with h5py.File('Redshift_GX.h5') as f:
    y0 = f['a']
    z = np.array(y0)
print(z)
#read the main_tree text(which saved in the Merg_tree_read file)
with h5py.File('main_tree_GX.h5') as f:
    y = f['a']
    main_tree = np.array(y)
print(main_tree[0,:])