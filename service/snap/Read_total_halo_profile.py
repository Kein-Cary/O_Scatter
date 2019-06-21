"""
###this file try to Read out all of the mass evolution of main halo
##from the redshift of the merger start
"""
import matplotlib as mpl
mpl.use('Agg')
from pygadgetreader import *
import pygadgetreader as pygdr
import numpy as np
import h5py
import pandas as pd
import astropy.io.ascii as asc
#import scipy.stats as st
#import matplotlib.pyplot as plt
#from handy import scatter as hsc
import find
import changds
#import glob
#def tot_halo_profile(ip):
_id_ = 0
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)
#read the meger tree to get the star time of each halo
with h5py.File('/home/cxkttwl/Scatter/MUSIC/Redshift.h5') as f:
    y0 = f['a']
    com_z = np.array(y0)
with h5py.File('/home/cxkttwl/Scatter/MUSIC/main_tree.h5') as f:
    y1 = f['a']
    tree_line = np.array(y1)
with h5py.File('/home/cxkttwl/Scatter/GadgetX/Redshift_GX.h5') as f:
    y2 = f['a']
    com_z_gx = np.array(y2)
with h5py.File('/home/cxkttwl/Scatter/GadgetX/main_tree_GX.h5') as f:
    y3 = f['a']
    tree_line_gx = np.array(y3)
L = tree_line.shape[0]
iv = np.zeros(L,dtype = np.int0)
for l in range(L):
    u = find.find1d(tree_line[l,:],0)
    iv[l] = u
L_gx = tree_line_gx.shape[0]
iv_gx = np.zeros(L_gx,dtype = np.int0)
for ll in range(L_gx):
    u_gx = find.find1d(tree_line_gx[ll,:],0)
    iv_gx[ll] = u_gx
id_z = com_z[iv[_id_]]
id_z_gx = com_z_gx[iv_gx[_id_]]
#read the data of MUSIC simulation

for k in range(len(alpha)):
    dd = alpha[-k]
    ss = str(dd)
    if len(ss) == 1:
        s_u = str('00%.0f' % dd)
    elif len(ss) == 2:
        s_u = str('0%.0f' % dd)
    else:
        s_u = str(dd)
    No_snap = s_u
    snap_z = pygdr.readheader(
            '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
    snap_z = np.abs(snap_z)
    print('Now redshift is %.3f' % snap_z) 
    if snap_z <= id_z:
        main_halo = asc.read(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z,
                converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
        goal_z = changds.inv_chidas('%.3f'%snap_z)
        goalz = find.find1d(com_z,goal_z)
        goal_halo = tree_line[_id_,goalz]
        #check the halo in this redshift
        try:
            HALO = np.array(main_halo['col1'])
            check_id = find.find1d(HALO,goal_halo) 
            Nbins = np.array(main_halo['col37'])
            # read the profile data
            halo_profile = pd.read_table(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%snap_z,
                    dtype = np.float)
            r = np.array(halo_profile['#r(1)'])
            r = np.abs(r)#to make sure that r is positive
            n_part = np.array(halo_profile['npart(2)'])
            m_in_r = np.array(halo_profile['M_in_r(3)'])
            m_star = np.array(halo_profile['M_star(26)'])
            m_gas = np.array(halo_profile['M_gas(25)'])
            dens = np.array(halo_profile['dens(5)'])
            ovdens = np.array(halo_profile['ovdens(4)'])
            local = np.int(Nbins[check_id]) 
            if check_id == 0:
                mg = m_gas[0:local]
                ms = m_star[0:local]
                mb = mg + ms
                md = m_in_r[0:local] - mb
                R = r[0:local]
            else:
                sum_bin = np.sum(Nbins[0:check_id])
                sum_bin = np.int(sum_bin)
                mg = m_gas[sum_bin:sum_bin+local]
                ms = m_star[sum_bin:sum_bin+local]
                mb = mg + ms
                md = m_in_r[sum_bin:sum_bin+local] - mb
                R = r[sum_bin:sum_bin+local]
            m_MU = np.array([mg, ms, mb, md, R])
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_halo-%.0f-profile_MU-%.3f.h5'%(k,goal_z),'w') as f:
                f['a'] = np.array(m_MU)
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_halo-%.0f-profile_MU-%.3f.h5'%(k,goal_z)) as f:
                for t in range(m_MU.shape[0]):
                    f['a'][t,:] = m_MU[t,:]  
        except IndexError:
            print('Now there is no merger tree to trace back!')
            break
for k in range(len(alpha)):
    dd = alpha[-k]
    ss = str(dd)
    if len(ss) == 1:
        s_u = str('00%.0f' % dd)
    elif len(ss) == 2:
        s_u = str('0%.0f' % dd)
    else:
        s_u = str(dd)
    No_snap = s_u
    snap_z_gx = pygdr.readheader(
            '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
    snap_z_gx = np.abs(snap_z_gx)
    print('Now redshift is %.3f' % snap_z_gx)
    if snap_z_gx <= id_z_gx:
        main_halo_gx = asc.read(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z_gx,
                converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
        goal_z_gx = changds.inv_chidas('%.3f'%snap_z_gx)
        goalz_gx = find.find1d(com_z_gx,goal_z_gx)
        goal_halo_gx = tree_line_gx[_id_,goalz_gx]
        #check the halo in this redshift
        try:
            HALO_gx = np.array(main_halo_gx['col1'])
            check_id_gx = find.find1d(HALO_gx,goal_halo_gx) 
            Nbins_gx = np.array(main_halo_gx['col37'])
            # read the profile data
            halo_profile_gx = pd.read_table(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%snap_z_gx,
                    dtype = np.float)
            r_gx = np.array(halo_profile_gx['#r(1)'])
            r_gx = np.abs(r_gx)#to make sure that r is positive
            n_part_gx = np.array(halo_profile_gx['npart(2)'])
            m_in_r_gx = np.array(halo_profile_gx['M_in_r(3)'])
            m_star_gx = np.array(halo_profile_gx['M_star(26)'])
            m_gas_gx = np.array(halo_profile_gx['M_gas(25)'])
            dens_gx = np.array(halo_profile_gx['dens(5)'])
            ovdens_gx = np.array(halo_profile_gx['ovdens(4)'])
            local_gx = np.int(Nbins_gx[check_id_gx]) 
            if check_id_gx == 0:
                mg_gx = m_gas_gx[0:local_gx]
                ms_gx = m_star_gx[0:local_gx]
                mb_gx = mg_gx + ms_gx
                md_gx = m_in_r_gx[0:local_gx] - mb_gx
                R_gx = r_gx[0:local_gx]
            else:
                sum_bin_gx = np.sum(Nbins_gx[0:check_id_gx])
                sum_bin_gx = np.int(sum_bin_gx)
                mg_gx = m_gas_gx[sum_bin_gx:sum_bin_gx+local_gx]
                ms_gx = m_star_gx[sum_bin_gx:sum_bin_gx+local_gx]
                mb_gx = mg_gx + ms_gx
                md_gx = m_in_r_gx[sum_bin_gx:sum_bin_gx+local_gx] - mb_gx  
                R_gx = r_gx[sum_bin_gx:sum_bin_gx+local_gx]
            m_GX = np.array([mg_gx, ms_gx, mb_gx, md_gx, R_gx])
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_halo-%.0f-profile_GX-%.3f.h5'%(k,goal_z_gx),'w') as f:
                f['a'] = np.array(m_GX)
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_halo-%.0f-profile_GX-%.3f.h5'%(k,goal_z_gx)) as f:
                for t in range(m_GX.shape[0]):
                    f['a'][t,:] = m_GX[t,:]    
        except IndexError:
            print('Now there is no merger tree to trace back!')
            break
#    returns
#tot_halo_profile(ip == True)