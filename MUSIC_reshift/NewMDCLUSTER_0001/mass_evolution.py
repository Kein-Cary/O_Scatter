import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
####part 0:read the redshift,to show the result in different z
with h5py.File('Redshift.h5') as f:
    y0 = f['a']
    z = np.array(y0)
with h5py.File('main_tree.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)
#get those halo which in range z<=1.0
ia = z<=1.0 #set the goal range of redshift
redshift = np.array(z[ia])
l = len(redshift)
halo_list = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
halo_id = np.array(halo_list['col1'])
_id_ = 0
goal_id = halo_id[_id_]#set the goal halo at z==0
#set the array for mass
halo_mass = np.zeros(l, dtype = np.float)
star_mass = np.zeros(l, dtype = np.float)
gas_mass = np.zeros(l, dtype = np.float)
dark_m = np.zeros(l,dtype = np.float)
for k in range(l):
    halo_read = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%redshift[k],
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
    id_read = np.array(halo_read['col1'])
    m_h_read = np.array(halo_read['col4'])#mass of halo
    m_s_read = np.array(halo_read['col65'])#mass of star
    m_g_read = np.array(halo_read['col45'])#mass of gas
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
    star_mass[k] = m_s_read[loc_on]
    gas_mass[k] = m_g_read[loc_on]
dark_m = halo_mass - star_mass - gas_mass
baryon_m = star_mass + gas_mass
plt.plot(redshift,halo_mass,'m-',label = r'$M_h$')
plt.plot(redshift,star_mass,'r-',label = r'$M_{\ast}$')
plt.plot(redshift,gas_mass,'g-',label = r'$M_g$')
plt.plot(redshift,baryon_m,'b-',label = r'$M_b$')
plt.plot(redshift,dark_m,'k-',label = r'$M_d$')
plt.legend(loc = 1)
plt.xlabel(r'$ z $')
plt.ylabel(r'$mass [M_\odot /h]$')
plt.yscale('log')
plt.title('GadgetMUSIC_mass')
#plt.savefig('mass_evolution_MUSIC',dpi=600)
plt.show()