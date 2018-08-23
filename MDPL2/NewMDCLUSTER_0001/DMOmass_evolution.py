import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
#'D:/mask/MDPL2/NewMDCLUSTER_0001/'
####part 0:read the redshift,to show the result in different z
with h5py.File('Redshift_DMO.h5') as f:
    y0 = f['a']
    z = np.array(y0)
with h5py.File('main_tree_DMO.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)
#get those halo which in range z<=1.0
ia = z<=3.0 #set the goal range of redshift
redshift = np.array(z[ia])
l = len(redshift)
halo_list = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z0.000.AHF_halos')
halo_id = np.array(halo_list['#ID(1)'])
_id_ = 0
goal_id = halo_id[_id_]#set the goal halo at z==0
#set the array for mass
halo_mass = np.zeros(l, dtype = np.float)
star_mass = np.zeros(l, dtype = np.float)
gas_mass = np.zeros(l, dtype = np.float)
for k in range(l):
    halo_read = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z%.3f.AHF_halos'%redshift[k])
    id_read = np.array(halo_read['#ID(1)'])
    m_h_read = np.array(halo_read['Mvir(4)'])#mass of halo
    ##identical the goal halo 
    ix = z == redshift[k]
    iy = ix.tolist()
    ID = iy.index(True)
    goal_halo = main_tree[_id_,ID]#get the goal halo at z!=0
    #get the properties
    ia = id_read == goal_halo
    ib = ia.tolist()
    loc_on = ib.index(True)
    halo_mass[k] = m_h_read[loc_on]
plt.plot(redshift,halo_mass,'b-',label = 'halo mass')
plt.legend(loc = 1)
plt.xlabel(r'$ z $')
plt.ylabel(r'$mass [M_\odot /h]$')
plt.yscale('log')
plt.title('DMO_mass')
plt.savefig('mass_evolution_DMO',dpi=600)
plt.show()