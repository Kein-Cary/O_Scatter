import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
import matplotlib.pyplot as plt
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
#figure the mass evolutionss
#def ma_evolution(idx):
idx = 0
with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/redshift.h5') as f:
    y0 = f['a']
    z = np.array(y0)
with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/main_tree.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)
#get those halo which in range z<=1.0
ia = z<=1.0 #set the goal range of redshift
redshift = np.array(z[ia])
l = len(redshift)
halo_list = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
halo_id = np.array(halo_list['col1'])
_id_ = idx
goal_id = halo_id[_id_]#set the goal halo at z==0
#set the array for mass
halo_mass = np.zeros(l, dtype = np.float)
star_mass = np.zeros(l, dtype = np.float)
gas_mass = np.zeros(l, dtype = np.float)
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
    iA = id_read == goal_halo
    iB = iA.tolist()
    loc_on = iB.index(True)
    halo_mass[k] = m_h_read[loc_on]
    star_mass[k] = m_s_read[loc_on]
    gas_mass[k] = m_g_read[loc_on]
baryon_m = star_mass + gas_mass
dark_mass = halo_mass - baryon_m 
mass_MUSIC = np.array([halo_mass,star_mass,gas_mass,baryon_m,dark_mass])
M = mass_MUSIC.shape[0]
with h5py.File('total_mass_MUSIC.h5','w') as f:
     f['a'] = np.array(mass_MUSIC)
with h5py.File('total_mass_MUSIC.h5') as f:  
    for t in range(M):
        f['a'][t,:] = mass_MUSIC[t,:]
###get the Data of GX
f = h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/Redshift_GX.h5','r')
y2 = f['a']
z_gx = np.array(y2)
f.close()
f = h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/main_tree_GX.h5','r') 
y3 = f['a']
main_tree_gx = np.array(y3)
f.close()
#get those halo which in range z<=1.0
ia_gx = z_gx<=1.0 #set the goal range of redshift
redshift_gx = np.array(z_gx[ia_gx])
l_gx = len(redshift_gx)
halo_list_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
halo_id_gx = np.array(halo_list_gx['col1'])
_id_gx = idx
goal_id_gx = halo_id_gx[_id_gx]#set the goal halo at z==0
#set the array for mass
halo_mass_gx = np.zeros(l_gx, dtype = np.float)
star_mass_gx = np.zeros(l_gx, dtype = np.float)
gas_mass_gx = np.zeros(l_gx, dtype = np.float)
for t in range(l_gx):
    halo_read_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%redshift_gx[t],
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
    id_read_gx = np.array(halo_read_gx['col1'])
    m_h_read_gx = np.array(halo_read_gx['col4'])#mass of halo
    m_s_read_gx = np.array(halo_read_gx['col65'])#mass of star
    m_g_read_gx = np.array(halo_read_gx['col45'])#mass of gas
    ##identical the goal halo 
    ix2 = z_gx == redshift_gx[t]
    iy2 = ix2.tolist()
    ID2 = iy2.index(True)
    goal_halo_gx = main_tree_gx[_id_gx,ID2]#get the goal halo at z!=0
    #get the properties
    iA2 = id_read_gx == goal_halo_gx
    iB2 = iA2.tolist()
    loc_on_gx = iB2.index(True)
    halo_mass_gx[t] = m_h_read_gx[loc_on_gx]
    star_mass_gx[t] = m_s_read_gx[loc_on_gx]
    gas_mass_gx[t] = m_g_read_gx[loc_on_gx]
baryon_m_gx = star_mass_gx + gas_mass_gx
dark_mass_gx = halo_mass_gx - baryon_m_gx
mass_GX = np.array([halo_mass_gx,star_mass_gx,gas_mass_gx,baryon_m_gx,dark_mass_gx])
N = mass_GX.shape[0]
with h5py.File('total_mass_GX.h5','w') as f:
     f['a'] = np.array(mass_GX)
with h5py.File('total_mass_GX.h5') as f:  
    for p in range(N):
        f['a'][p,:] = mass_GX[p,:]
plt.figure(figsize=(8,9))
plt.plot(redshift,halo_mass,'b-',label = r'$M_h[MUSIC]$')
plt.plot(redshift,star_mass,'r-',label = r'$M_s[MUSIC]$')
plt.plot(redshift,gas_mass,'g-',label = r'$M_g[MUSIC]$')
plt.plot(redshift,baryon_m,'m-',label = r'$M_b[MUSIC]$')
plt.plot(redshift,dark_mass,'k-',label = r'$M_d[MUSIC]$')
plt.plot(redshift_gx,halo_mass_gx,'b--',label = r'$M_h[GX]$')
plt.plot(redshift_gx,star_mass_gx,'r--',label = r'$M_s[GX]$')
plt.plot(redshift_gx,gas_mass_gx,'g--',label = r'$M_g[GX]$')
plt.plot(redshift_gx,baryon_m_gx,'m--',label = r'$M_b[GX]$')
plt.plot(redshift_gx,dark_mass_gx,'k--',label = r'$M_d[GX]$')
plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=2,
   ncol=5, mode = "expand",fontsize=7.5, borderaxespad=1.) 
plt.xlabel(r'$ z $')
plt.ylabel(r'$M [M_\odot /h]$')
plt.yscale('log')
plt.title('M-z relation')
plt.savefig('M-z',dpi=600)
plt.show()

f = plt.figure(figsize = (8,9))
f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
ratio_gas_b = gas_mass/baryon_m
ratio_gas_d = gas_mass/dark_mass
ratio_star_d = star_mass/dark_mass
ratio_star_b = star_mass/baryon_m
ratio_b = baryon_m/halo_mass
ratio_d = dark_mass/halo_mass
ratio_gas_b_gx = gas_mass_gx/baryon_m_gx
ratio_gas_d_gx = gas_mass_gx/dark_mass_gx
ratio_star_d_gx = star_mass_gx/dark_mass_gx
ratio_star_b_gx = star_mass_gx/baryon_m_gx
ratio_b_gx = baryon_m_gx/halo_mass_gx
ratio_d_gx = dark_mass_gx/halo_mass_gx

ax1.plot(redshift,ratio_gas_b,'b-',label = r'$\eta -M_g-M_b$')
ax1.plot(redshift,ratio_d,'k-',label = r'$\eta -M_d-M_h$')
ax1.plot(redshift_gx,ratio_gas_b_gx,'b--',label = r'$\eta -M_g-M_b -GX$')
ax1.plot(redshift_gx,ratio_d_gx,'k--',label = r'$\eta -M_d-M_h -GX$')
ax1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=2,
     ncol=2, fontsize=7.5, borderaxespad=0.) 
plt.ylabel(r'$\eta$')
plt.ylim(0.75,0.95)
ax1.set_title(r'$\eta - z $')

ax2.plot(redshift,ratio_gas_d,'g-',label = r'$\eta[M_g-M_d]$')
ax2.plot(redshift,ratio_star_d,'r-',label = r'$\eta[M_\ast-M_d]$')
ax2.plot(redshift,ratio_star_b,'m-',label = r'$\eta[M_\ast-M_b]$')
ax2.plot(redshift,ratio_b,'y-',label = r'$\eta[M_b-M_h]$')
ax2.plot(redshift_gx,ratio_gas_d_gx,'g--',label = r'$\eta[M_g-M_d]GX$')
ax2.plot(redshift_gx,ratio_star_d_gx,'r--',label = r'$\eta[M_\ast-M_d]GX$')
ax2.plot(redshift_gx,ratio_star_b_gx,'m--',label = r'$\eta[M_\ast-M_b]GX$')
ax2.plot(redshift_gx,ratio_b_gx,'y--',label = r'$\eta[M_b-M_h]GX$')
ax2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=2,
     ncol=4, fontsize=7.5, borderaxespad=0.) 
plt.ylim(0.,0.3)
plt.ylabel(r'$\eta$')
plt.xlabel(r'$ z $')
#plt.title(r'$\eta-z relation$')
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.savefig('ratio of mass',dpi=600)
plt.show()