#This file use for saving these data read from the simulation
#for I/O quakly,use H5py files to save
import sys
sys.path.insert(0,'D:/mask')
#import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import glob
import h5py 
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
#import matplotlib.pyplot as plt
##firstly,try to read the red_shift value

namestr = glob.glob('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001*AHF_halos')
#print(namestr[0])
red_shift = [float(namestr[70:-10]) for namestr in namestr]
print(red_shift)

##sort to red_shift and save as h5py file
redshift = np.sort(red_shift)
with h5py.File('redshift_GX.h5','w') as f:
     f['a'] = np.array(redshift)
#print(redshift)
###just one time to save the redshift data

with h5py.File('redshift_GX.h5') as f:
     redshift = np.array(f['a'])
##next step,try to read the analysis catalogue from the data,and save as h5py files
###firstly,read the halo mass and the red_shift
##get the number of main halo
main_value = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_mtree_idx',
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
N = len(main_value['col1'])
M = len(redshift)
main_tree = np.zeros((N,M),dtype = np.int64)#save the link to main progenitor
with h5py.File('main_tree_GX.h5','w') as f:
     f['a'] = np.array(main_tree)
for k in range(0,len(redshift)-1):
    if redshift[k] <= 9.85 :
        id_value = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_mtree_idx'%redshift[k],
                    converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
        if k == 0 :
            main_tree[:,0] = id_value['col1']
            main_tree[:,1] = id_value['col2']
        else:
            for t in range(len(id_value['col1'])):
                ip_halo = main_tree[:,k] == id_value['col1'][t]
                ipx = ip_halo.tolist()
                ipy = True in ipx
                if ipy == True:
                    T_ip = ipx.index(True)
                    main_tree[T_ip,k+1] = id_value['col2'][t]
with h5py.File('main_tree_GX.h5') as f:
     for k in range(len(main_tree)):
         f['a'][k,:] = main_tree[k,:]
''' 
#next file try to think about all the halo(involved sub-halo)           
basic_tree = np.zeros((M,N),dtype = np.int64)#save all the evolution    
with h5py.File('basic_tree_GX.h5','w') as f:
     f['a'] = np.array(basic_tree)
'''