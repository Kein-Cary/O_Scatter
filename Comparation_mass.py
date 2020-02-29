import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
import matplotlib.pyplot as plt

idx = 0
with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/redshift.h5') as f:
    y0 = f['a']
    z = np.array(y0)
with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/main_tree.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)

l = len(z)
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
    if ((z[k] >= 3.20) & (z[k] <= 3.97)) | ((z[k] >= 4.30) & (z[k] <= 5.00)):
        halo_read = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%z[k],
                    converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
        id_read = np.array(halo_read['col1'])
        sub_id = np.array(halo_read['col2'])
        m_h_read = np.array(halo_read['col4'])#mass of halo
        m_s_read = np.array(halo_read['col65'])#mass of star
        m_g_read = np.array(halo_read['col45'])#mass of gas

        Rvir = np.array(halo_read['col12'])
        xhalo = np.array(halo_read['col6'])
        yhalo = np.array(halo_read['col7'])
        zhalo = np.array(halo_read['col8'])
        raise
        ##identical the goal halo 
        ix = z == z[k]
        iy = ix.tolist()
        ID = iy.index(True)
        goal_halo = main_tree[_id_,ID]#get the goal halo at z!=0
        #get the properties
        iA = id_read == goal_halo
        iB = iA.tolist()
        loc_on = iB.index(True)

        print('k =', k)
        print('*'*20)
        halo_mass[k] = m_h_read[loc_on]
        star_mass[k] = m_s_read[loc_on]
        gas_mass[k] = m_g_read[loc_on]

        x0 = xhalo[loc_on]
        y0 = yhalo[loc_on]
        z0 = zhalo[loc_on]
        R0 = Rvir[loc_on]

        print('DM = ', halo_mass[k])
        print('star = ', star_mass[k])
        print('gas = ', gas_mass[k])
        # satellite choose
        ih = sub_id == id_read[loc_on]
        msat = m_s_read[ih]
        print('sat_total = ', msat)
        cen_ip = np.where(msat == np.max(msat))[0]
        c_id = np.where(ih == True)[0]
        cen_mass = m_s_read[c_id[cen_ip[0]]]
        print('m_cen = ', cen_mass)
        cen_x = xhalo[c_id[cen_ip[0]]]
        cen_y = yhalo[c_id[cen_ip[0]]]
        cen_z = zhalo[c_id[cen_ip[0]]]

        sub_x = xhalo[ih]
        sub_y = yhalo[ih]
        sub_z = zhalo[ih]
        sub_r = Rvir[ih]

        dr_sub = np.sqrt((sub_x - x0)**2 
            + (sub_y - y0)**2 + (sub_z - z0)**2)
        isub = dr_sub <= R0

        real_sub = msat[isub]
        real_x = sub_x[isub]
        real_y = sub_y[isub]
        real_z = sub_z[isub]

        dr_sub = np.sqrt((sub_x - x0)**2 
            + (sub_y - y0)**2 + (sub_z - z0)**2)
        isub = dr_sub <= R0
        real_sub = msat[isub]
        print('sat_in_200 =', real_sub)
        d_msat = np.sum(real_sub) - cen_mass
        print('msat = ', d_msat)
        print('*'*20)

        plt.figure( figsize = (16, 5))
        ax1 = plt.subplot(131)
        ax2 = plt.subplot(132)
        ax3 = plt.subplot(133)

        ax1.plot(x0, y0, '+', c = 'b')
        ax1.scatter(sub_x, sub_y, s = 15, marker = 'o', facecolors = '', edgecolors = 'g', alpha = 0.5, linewidths = 1)
        ax1.scatter(real_x, real_y, s = 5, marker = 'o', facecolors = 'r',  edgecolors = '', alpha = 0.5, linewidths = 1)
        ax1.scatter(cen_x, cen_y, s = 10, marker = 'X', c = 'k', linewidths = 1)
        disk1 = plt.Circle((x0, y0), R0, color='r', fill = False)
        ax1.add_artist(disk1)

        ax1.set_title('x-y')
        ax1.set_xlabel('r[kpc/h]')
        ax1.set_ylabel('r[kpc/h]')
        ax1.axis('equal')
        ax1.set_xlim(x0-R0, x0+R0)
        ax1.set_ylim(y0-R0, y0+R0)
        ax1.tick_params(axis = 'both', which = 'both', direction = 'in')

        ax2.plot(x0, z0, '+', c = 'b')
        ax2.scatter(sub_x, sub_z, s = 15, marker = 'o',  facecolors = '', edgecolors = 'g', alpha = 0.5, linewidths = 1)
        ax2.scatter(real_x, real_z, s = 5, marker = 'o', facecolors = 'r',  edgecolors = '', alpha = 0.5, linewidths = 1)
        ax2.scatter(cen_x, cen_z, s = 10, marker = 'X', c = 'k', linewidths = 1)
        disk2 = plt.Circle((x0, z0), R0, color='r', fill = False)
        ax2.add_artist(disk2)

        ax2.set_title('x-z')
        ax2.set_xlabel('r[kpc/h]')
        ax2.set_ylabel('r[kpc/h]')
        ax2.axis('equal')
        ax2.set_xlim(x0-R0, x0+R0)
        ax2.set_ylim(z0-R0, z0+R0)
        ax2.tick_params(axis = 'both', which = 'both', direction = 'in')

        ax3.plot(y0, z0, '+', c = 'b')
        ax3.scatter(sub_y, sub_z, s = 15, marker = 'o', facecolors = '', edgecolors = 'g', alpha = 0.5, linewidths = 1)
        ax3.scatter(real_y, real_z, s = 5, marker = 'o', facecolors = 'r',  edgecolors = '', alpha = 0.5, linewidths = 1)
        ax3.scatter(cen_y, cen_z, s = 10, marker = 'X', c = 'k', linewidths = 1)
        disk3 = plt.Circle((y0, z0), R0, color='r', fill = False)
        ax3.add_artist(disk3)

        ax3.set_title('y-z')
        ax3.set_xlabel('r[kpc/h]')
        ax3.set_ylabel('r[kpc/h]')
        ax3.axis('equal')
        ax3.set_xlim(y0-R0, y0+R0)
        ax3.set_ylim(z0-R0, z0+R0)
        ax3.tick_params(axis = 'both', which = 'both', direction = 'in')

        plt.tight_layout()
        plt.savefig('MU_test_show_z%.3f.png'%z[k], dpi = 300)
        plt.close()
'''
baryon_m = star_mass + gas_mass
dark_mass = halo_mass - baryon_m 
mass_MUSIC = np.array([halo_mass,star_mass,gas_mass,baryon_m,dark_mass])
M = mass_MUSIC.shape[0]
with h5py.File('total_mass_MUSIC.h5','w') as f:
     f['a'] = np.array(mass_MUSIC)
with h5py.File('total_mass_MUSIC.h5') as f:  
    for t in range(M):
        f['a'][t,:] = mass_MUSIC[t,:]
'''
###get the Data of GX
f = h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/Redshift_GX.h5','r')
y2 = f['a']
z_gx = np.array(y2)
f.close()
f = h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/main_tree_GX.h5','r') 
y3 = f['a']
main_tree_gx = np.array(y3)
f.close()

l_gx = len(z_gx)
halo_list_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
halo_id_gx = np.array(halo_list_gx['col1'])
_id_gx = idx
goal_id_gx = halo_id_gx[_id_gx]#set the goal halo at z==0
#set the array for mass
halo_mass = np.zeros(l_gx, dtype = np.float)
star_mass = np.zeros(l_gx, dtype = np.float)
gas_mass = np.zeros(l_gx, dtype = np.float)
for k in range(l_gx):
    if (z_gx[k] >= 2.60) & (z_gx[k] <= 3.20):
        halo_read_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'% z_gx[k],
                    converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
        id_read = np.array(halo_read_gx['col1'])
        sub_id = np.array(halo_read_gx['col2'])
        m_h_read = np.array(halo_read_gx['col4'])#mass of halo
        m_s_read = np.array(halo_read_gx['col65'])#mass of star
        m_g_read = np.array(halo_read_gx['col45'])#mass of gas

        Rvir = np.array(halo_read_gx['col12'])
        xhalo = np.array(halo_read_gx['col6'])
        yhalo = np.array(halo_read_gx['col7'])
        zhalo = np.array(halo_read_gx['col8'])

        ##identical the goal halo 
        ix2 = z_gx == z_gx[k]
        iy2 = ix2.tolist()
        ID2 = iy2.index(True)
        goal_halo = main_tree_gx[_id_gx,ID2]#get the goal halo at z!=0
        #get the properties
        iA2 = id_read == goal_halo
        iB2 = iA2.tolist()
        loc_on = iB2.index(True)

        print('k =', k)
        print('*'*20)
        halo_mass[k] = m_h_read[loc_on]
        star_mass[k] = m_s_read[loc_on]
        gas_mass[k] = m_g_read[loc_on]

        x0 = xhalo[loc_on]
        y0 = yhalo[loc_on]
        z0 = zhalo[loc_on]
        R0 = Rvir[loc_on]

        print('DM = ', halo_mass[k])
        print('star = ', star_mass[k])
        print('gas = ', gas_mass[k])
        # satellite choose
        ih = sub_id == id_read[loc_on]
        msat = m_s_read[ih]
        print('sat_total = ', msat)
        cen_ip = np.where(msat == np.max(msat))[0]
        c_id = np.where(ih == True)[0]
        cen_mass = m_s_read[c_id[cen_ip[0]]]
        print('m_cen = ', cen_mass)
        cen_x = xhalo[c_id[cen_ip[0]]]
        cen_y = yhalo[c_id[cen_ip[0]]]
        cen_z = zhalo[c_id[cen_ip[0]]]

        sub_x = xhalo[ih]
        sub_y = yhalo[ih]
        sub_z = zhalo[ih]
        sub_r = Rvir[ih]

        dr_sub = np.sqrt((sub_x - x0)**2 
            + (sub_y - y0)**2 + (sub_z - z0)**2)
        isub = dr_sub <= R0
        real_sub = msat[isub]

        dr_sub = np.sqrt((sub_x - x0)**2 
            + (sub_y - y0)**2 + (sub_z - z0)**2)
        isub = dr_sub <= R0
        real_sub = msat[isub]
        real_x = sub_x[isub]
        real_y = sub_y[isub]
        real_z = sub_z[isub]
        print('sat_in_200 =', real_sub)
        d_msat = np.sum(real_sub) - cen_mass
        print('msat = ', d_msat)
        print('*'*20)

        plt.figure( figsize = (16, 5))
        ax1 = plt.subplot(131)
        ax2 = plt.subplot(132)
        ax3 = plt.subplot(133)

        ax1.plot(x0, y0, '+', c = 'b')
        ax1.scatter(sub_x, sub_y, s = 15, marker = 'o', facecolors = '', edgecolors = 'g', alpha = 0.5, linewidths = 1)
        ax1.scatter(real_x, real_y, s = 5, marker = 'o', facecolors = 'r',  edgecolors = '', alpha = 0.5, linewidths = 1)
        ax1.scatter(cen_x, cen_y, s = 10, marker = 'X', c = 'k', linewidths = 1)
        disk1 = plt.Circle((x0, y0), R0, color='r', fill = False)
        ax1.add_artist(disk1)

        ax1.set_title('x-y')
        ax1.set_xlabel('r[kpc/h]')
        ax1.set_ylabel('r[kpc/h]')
        ax1.axis('equal')
        ax1.set_xlim(x0-R0, x0+R0)
        ax1.set_ylim(y0-R0, y0+R0)
        ax1.tick_params(axis = 'both', which = 'both', direction = 'in')

        ax2.plot(x0, z0, '+', c = 'b')
        ax2.scatter(sub_x, sub_z, s = 15, marker = 'o',  facecolors = '', edgecolors = 'g', alpha = 0.5, linewidths = 1)
        ax2.scatter(real_x, real_z, s = 5, marker = 'o', facecolors = 'r',  edgecolors = '', alpha = 0.5, linewidths = 1)
        ax2.scatter(cen_x, cen_z, s = 10, marker = 'X', c = 'k', linewidths = 1)
        disk2 = plt.Circle((x0, z0), R0, color='r', fill = False)
        ax2.add_artist(disk2)

        ax2.set_title('x-z')
        ax2.set_xlabel('r[kpc/h]')
        ax2.set_ylabel('r[kpc/h]')
        ax2.axis('equal')
        ax2.set_xlim(x0-R0, x0+R0)
        ax2.set_ylim(z0-R0, z0+R0)
        ax2.tick_params(axis = 'both', which = 'both', direction = 'in')

        ax3.plot(y0, z0, '+', c = 'b')
        ax3.scatter(sub_y, sub_z, s = 15, marker = 'o', facecolors = '', edgecolors = 'g', alpha = 0.5, linewidths = 1)
        ax3.scatter(real_y, real_z, s = 5, marker = 'o', facecolors = 'r',  edgecolors = '', alpha = 0.5, linewidths = 1)
        ax3.scatter(cen_y, cen_z, s = 10, marker = 'X', c = 'k', linewidths = 1)
        disk3 = plt.Circle((y0, z0), R0, color='r', fill = False)
        ax3.add_artist(disk3)

        ax3.set_title('y-z')
        ax3.set_xlabel('r[kpc/h]')
        ax3.set_ylabel('r[kpc/h]')
        ax3.axis('equal')
        ax3.set_xlim(y0-R0, y0+R0)
        ax3.set_ylim(z0-R0, z0+R0)
        ax3.tick_params(axis = 'both', which = 'both', direction = 'in')

        plt.tight_layout()
        plt.savefig('GX_test_show_z%.3f.png'%z_gx[k], dpi = 300)
        plt.close()
raise
'''
baryon_m_gx = star_mass_gx + gas_mass_gx
dark_mass_gx = halo_mass_gx - baryon_m_gx
mass_GX = np.array([halo_mass_gx,star_mass_gx,gas_mass_gx,baryon_m_gx,dark_mass_gx])
N = mass_GX.shape[0]
with h5py.File('total_mass_GX.h5','w') as f:
     f['a'] = np.array(mass_GX)
with h5py.File('total_mass_GX.h5') as f:  
    for p in range(N):
        f['a'][p,:] = mass_GX[p,:]
'''