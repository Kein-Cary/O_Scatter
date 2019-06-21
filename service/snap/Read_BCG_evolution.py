# this file try to find the BCG in a given halo
import matplotlib as mpl
mpl.use('Agg')
from pygadgetreader import *
import pygadgetreader as pygdr
import numpy as np
import h5py
import pandas as pd
import astropy.io.ascii as asc
import scipy.stats as st
import matplotlib.pyplot as plt
from handy import scatter as hsc
import find
import changds
def BCG_evolution(ip,resolution,g_size):
    """
    ip: the id of the main halo at z == 0
    resolution: the resollution to find the centre of BCG
    g_size:assum the size of BCG in radius,in unit kpc/h
    all of these parameter is useng for MUSIC and GadgetX simulation
    """
    resolution = resolution
    size_BCG = g_size
    _id_ = ip
    Nr = 1./resolution
    # to make sure the bins density is the same
    N_size = 50.*(size_BCG/100.0)
    Nsize = np.int0(np.ceil(N_size))
    alpha = np.linspace(0, 128, 129)
    alpha = np.int0(alpha)
    T_scale = len(alpha)
    # set array to save the data of BCGs:m--mass,rho--density profile,m-rho--mean density profile
    m_s = np.zeros(T_scale,dtype = np.float)
    m_g = np.zeros(T_scale,dtype = np.float)
    m_bh = np.zeros(T_scale,dtype = np.float)
    m_s_gx = np.zeros(T_scale,dtype = np.float)
    m_g_gx = np.zeros(T_scale,dtype = np.float)
    m_bh_gx = np.zeros(T_scale,dtype = np.float)
    rho_s = np.zeros((T_scale,Nsize-1),dtype = np.float)
    rho_g = np.zeros((T_scale,Nsize-1),dtype = np.float)
    rho_s_gx = np.zeros((T_scale,Nsize-1),dtype = np.float)
    rho_g_gx = np.zeros((T_scale,Nsize-1),dtype = np.float)
    m_rho_s = np.zeros((T_scale,Nsize),dtype = np.float)
    m_rho_g = np.zeros((T_scale,Nsize),dtype = np.float)
    m_rho_s_gx = np.zeros((T_scale,Nsize),dtype = np.float)
    m_rho_g_gx = np.zeros((T_scale,Nsize),dtype = np.float)
    r_mu = np.logspace(-2, np.log10(size_BCG), Nsize)
    r_gx = np.logspace(-2, np.log10(size_BCG), Nsize)
    z_mu = np.zeros(T_scale,dtype = np.float)
    z_gx = np.zeros(T_scale,dtype = np.float)
    r = np.array([r_mu,r_gx])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_R_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(r)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_R_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(np.int0(r.shape[0])):
            f['a'][t,:] = np.array(r[t,:])
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
    # section2:read out the position data of MUSIC
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
        snap_z = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
        snap_z = np.abs(snap_z)
        print('Now redshift is %.3f' % snap_z)
        snap_N = pygdr.readheader(
            '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
        print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f' % (
            snap_N[0], snap_N[1], snap_N[2], snap_N[3], snap_N[4], snap_N[5]))
        #snap_name = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'header')
        # for type
        id_z = com_z[iv[_id_]]
        if snap_z <= id_z:
            z_mu[k] = snap_z
            snap_shot_gas = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_mass_gas = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
            snap_shot_DM = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
            try:
                snap_shot_star = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
                snap_mass_star = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
            except SystemExit:print('no star particles now')
            # for respective position
            snap_shot_bulge = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bulge')
            snap_shot_disk = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'disk')
            try:
                snap_shot_bndry = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bndry')
                snap_mass_BH = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','bndry')
            except SystemExit:
                print('no boundary particles now')
            # for density profile of gas(there only density about gas in the simulation data)
            # snap_dens = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'rho','gas')
            # the total mass distribution of snap_shot_128
            # next,try to find the particles belong to halo 128000000000001(in MUSIC simulation)
            main_halo = asc.read(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z,
                                 converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z = changds.inv_chidas('%.3f'%snap_z)
            goalz = find.find1d(com_z,goal_z)
            goal_halo = tree_line[_id_,goalz]
            #check the halo in this redshift
            HALO = np.array(main_halo['col1'])
            check_id = find.find1d(HALO,goal_halo) 
            Rvir = np.array(main_halo['col12'])
            xhalo = np.array(main_halo['col6'])
            yhalo = np.array(main_halo['col7'])
            zhalo = np.array(main_halo['col8'])
            x0 = xhalo[check_id]
            y0 = yhalo[check_id]
            z0 = zhalo[check_id]
            R0 = Rvir[check_id]
            dgas = np.sqrt((snap_shot_gas[:, 0]-x0)**2 +
                           (snap_shot_gas[:, 1]-y0)**2 + (snap_shot_gas[:, 2]-z0)**2)
            ig = dgas <= R0
            inlgas = snap_shot_gas[ig, :]
            inlmass_gas = snap_mass_gas[ig]
            
            dDM = np.sqrt((snap_shot_DM[:, 0]-x0)**2 +
                          (snap_shot_DM[:, 1]-y0)**2 + (snap_shot_DM[:, 2]-z0)**2)
            iD = dDM <= R0
            inlDM = snap_shot_DM[iD, :]
            ddisk = np.sqrt((snap_shot_disk[:, 0]-x0)**2 +
                            (snap_shot_disk[:, 1]-y0)**2 + (snap_shot_disk[:, 2]-z0)**2)
            idisk = ddisk <= R0
            inldisk = snap_shot_disk[idisk, :]
            dbulge = np.sqrt((snap_shot_bulge[:, 0]-x0)**2 +
                             (snap_shot_bulge[:, 1]-y0)**2 + (snap_shot_bulge[:, 2]-z0)**2)
            ibu = dbulge <= R0
            inlbulge = snap_shot_bulge[ibu, :]
            
            if snap_N[-2] != 0:
                dstar = np.sqrt((snap_shot_star[:, 0]-x0)**2+(snap_shot_star[:, 1]-y0)**2 + (snap_shot_star[:, 2]-z0)**2)
                ids = dstar <= R0
                inlstar = snap_shot_star[ids, :]
                inlmass_star = snap_mass_star[ids]
            if snap_N[-1] != 0:
                dbndry = np.sqrt(
                    (snap_shot_bndry[:, 0]-x0)**2+(snap_shot_bndry[:, 1]-y0)**2 + (snap_shot_bndry[:, 2]-z0)**2)
                idb = dbndry <= R0
                inlbndry = snap_shot_bndry[idb, :]
                inlmass_BH = snap_mass_BH[idb]
            # try to find the most densitive area in the central area(which will set as BCG)
            R_range = 150.0
            # set the central area as a cubic
            edge_x = np.array([x0-R_range,x0+R_range])
            edge_y = np.array([y0-R_range,y0+R_range])
            edge_z = np.array([z0-R_range,z0+R_range])
            iv_s = ((inlstar[:,0]<=edge_x[1])&(inlstar[:,0]>=edge_x[0]))&(
                    (inlstar[:,1]<=edge_y[1])&(inlstar[:,1]>=edge_y[0]))&((inlstar[:,2]<=edge_z[1])&(inlstar[:,2]>=edge_z[0]))
            inl_star = inlstar[iv_s,:]
            # central postion and density of star
            num_bins = np.ceil(R_range*2/resolution)
            try:
                hist_star, edge_star = np.histogramdd(inl_star, bins=(num_bins, num_bins, num_bins))
            except ValueError:
                #print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
                break
            bin_x_star = np.array(edge_star[0])
            bin_y_star = np.array(edge_star[1])
            bin_z_star = np.array(edge_star[2])
            inumber1 = hist_star >= 10
            #to find the first five density center and then choose the closed one to the halo center
            maxN = len(inumber1)
            cen_po_star = np.zeros((maxN,3),dtype = np.float)
            hist_use = hist_star
            for p in range(maxN):
                is_max = np.unravel_index(np.argmax(hist_use, axis=None), hist_use.shape)
                cenxstar = (bin_x_star[is_max[0]+1]+bin_x_star[is_max[0]])/2.0
                cenystar = (bin_y_star[is_max[1]+1]+bin_y_star[is_max[1]])/2.0
                cenzstar = (bin_z_star[is_max[2]+1]+bin_z_star[is_max[2]])/2.0
                cen_po_star[p,:] = np.array([cenxstar,cenystar,cenzstar])
                hist_use[is_max] = 0.0
            compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+
                               (cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
            ismin = find.find1d(compare_d,np.min(compare_d))
            cen_x_star = cen_po_star[ismin,0]
            cen_y_star = cen_po_star[ismin,1]
            cen_z_star = cen_po_star[ismin,2]
            # next,give radius and calculate the properties,in units kpc/h
            """
            according to Milky Way,assuming the BCG radiu is about 30kpc/h,but for finding the cut point
            give a 100kpc/h,to find the mean density profile
            """
            R_BCG = size_BCG
            r_bcg = np.logspace(1e-3,np.log10(R_BCG),Nsize)
            n_r = len(r_bcg)
            inl_m_tot = np.zeros(n_r,dtype = np.float)
            inl_m_g = np.zeros(n_r,dtype = np.float)
            #above three mass instand:mass of a gas particle,star particle and DM particle aspectlly
            r_gas = np.sqrt((inlgas[:,0]-cen_x_star)**2+(inlgas[:,1]-cen_y_star)**2+
                            (inlgas[:,2]-cen_z_star)**2)
            for p in range(n_r):
                ia = r_gas <= r_bcg[p]
                inl_m_g[p] = np.sum(inlmass_gas[ia]*10**10)
            m_g[k] = inl_m_g[-1]
            if snap_N[-2] !=0:
                inl_m_s = np.zeros(n_r,dtype = np.float)
                r_star = np.sqrt((inlstar[:,0]-cen_x_star)**2+(inlstar[:,1]-cen_y_star)**2+
                             (inlstar[:,2]-cen_z_star)**2)
                for p in range(n_r):
                    ib = r_star <= r_bcg[p]
                    inl_m_s[p] = np.sum(inlmass_star[ib]*10**10)
                m_s[k] = inl_m_s[-1]
            if snap_N[-1] !=0:
                inl_BH = {}
                inl_BH_m = np.zeros(n_r,dtype = np.float)
                r_BH = np.sqrt((inlbndry[:,0]-cen_x_star)**2+(inlbndry[:,1]-cen_y_star)**2+
                             (inlbndry[:,2]-cen_z_star)**2)
                for p in range(n_r):
                    ic = r_BH <= r_bcg[p]
                    inl_BH[p] = inlbndry[ic,:]
                    inl_BH_m[p] = np.sum(inlmass_BH[ic]*10**10)
                m_bh[k] = inl_BH_m[-1]
            if snap_N[-1] !=0 and snap_N[-2] !=0:
                inl_m_tot = inl_m_g+ inl_m_s + inl_BH_m
            if snap_N[-1] ==0 and snap_N[-2] !=0:
                inl_m_tot = inl_m_g+ inl_m_s
            if snap_N[-1]==0 and snap_N[-2] ==0:
                inl_m_tot = inl_m_g
            inl_m_b = inl_m_g + inl_m_s
            inl_rho_b = np.zeros(n_r-1,dtype = np.float)
            inl_rho_g = np.zeros(n_r-1,dtype = np.float)
            inl_rho_s = np.zeros(n_r-1,dtype = np.float)
            for p in range(n_r-1):
                dr = r_bcg[p+1] - r_bcg[p]
                inl_rho_g[p] = (inl_m_g[p+1]-inl_m_g[p])/(4*np.pi*dr*r_bcg[p]**2)
                if snap_N[-2] !=0:
                    inl_rho_s[p] = (inl_m_s[p+1]-inl_m_s[p])/(4*np.pi*dr*r_bcg[p]**2)
            rho_g[k,:] = inl_rho_g
            rho_s[k,:] = inl_rho_s
            inl_rho_b = inl_rho_g + inl_rho_s
            mean_rho_g = inl_m_g/(4*np.pi*r_bcg**3/3)
            m_rho_g[k,:] = mean_rho_g
            if snap_N[-2] !=0:
                mean_rho_s = inl_m_s/(4*np.pi*r_bcg**3/3)
                m_rho_s[k,:] = mean_rho_s
            mean_rho_b = inl_m_b/(4*np.pi*r_bcg**3/3)
    # section 3: read the data of GX simution
    id_z_gx = com_z_gx[iv_gx[_id_]]
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
        snap_z_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
        snap_z_gx = np.abs(snap_z_gx)
        print('Now redshift is %.3f' % snap_z_gx)
        snap_N_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
        print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
    	      % (snap_N_gx[0], snap_N_gx[1], snap_N_gx[2], snap_N_gx[3], snap_N_gx[4], snap_N_gx[5]))
        #snap_name_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'header')
        # for type
        if snap_z_gx <= id_z_gx:
            z_gx[k] = snap_z_gx
            snap_shot_gas_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_mass_gas_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
            snap_shot_DM_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
            try:
                snap_shot_star_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
                snap_mass_star_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
            except SystemExit:
                print('no star particles now')
            # for respective position
            snap_shot_bulge_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bulge')
            snap_shot_disk_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'disk')
            try:
                snap_shot_bndry_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bndry')
                snap_mass_BH_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','bndry')
            except SystemExit:
                print('no boundary particles now')
            # for density profile of gas(there only density about gas in the simulation data)
            # snap_dens_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'rho','gas')
            # the total mass distribution of snap_shot_128
            # next,try to find the particles belong to halo 128000000000001(in GX simulation)
            main_halo_gx = asc.read(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z_gx,
                                converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z_gx = changds.inv_chidas('%.3f'%snap_z_gx)
            goalz_gx = find.find1d(com_z_gx,goal_z_gx)
            goal_halo_gx = tree_line_gx[_id_,goalz_gx]
            #check the halo in this redshift
            HALO_gx = np.array(main_halo_gx['col1'])
            check_id_gx = find.find1d(HALO_gx,goal_halo_gx) 
            Rvir_gx = np.array(main_halo_gx['col12'])
            xhalo_gx = np.array(main_halo_gx['col6'])
            yhalo_gx = np.array(main_halo_gx['col7'])
            zhalo_gx = np.array(main_halo_gx['col8'])
            x0_gx = xhalo_gx[check_id_gx]
            y0_gx = yhalo_gx[check_id_gx]
            z0_gx = zhalo_gx[check_id_gx]
            R0_gx = Rvir_gx[check_id_gx]
            dgas_gx = np.sqrt((snap_shot_gas_gx[:, 0]-x0_gx)**2+(snap_shot_gas_gx[:, 1]-y0_gx)**2 +(snap_shot_gas_gx[:, 2]-z0_gx)**2)
            ig_gx = dgas_gx <= R0_gx
            inlgas_gx = snap_shot_gas_gx[ig_gx, :]
            inlmass_gas_gx = snap_mass_gas_gx[ig_gx]
            
            dDM_gx = np.sqrt((snap_shot_DM_gx[:, 0]-x0_gx)**2+(snap_shot_DM_gx[:, 1]-y0_gx)**2 +(snap_shot_DM_gx[:, 2]-z0_gx)**2)
            iD_gx = dDM_gx <= R0_gx
            inlDM_gx = snap_shot_DM_gx[iD_gx, :]
            ddisk_gx = np.sqrt((snap_shot_disk_gx[:, 0]-x0_gx)**2+(snap_shot_disk_gx[:, 1]-y0_gx)**2 +(snap_shot_disk_gx[:, 2]-z0_gx)**2)
            idisk_gx = ddisk_gx <= R0_gx
            inldisk_gx = snap_shot_disk_gx[idisk_gx, :]
            dbulge_gx = np.sqrt((snap_shot_bulge_gx[:, 0]-x0_gx)**2+(snap_shot_bulge_gx[:, 1]-y0_gx)**2 +(snap_shot_bulge_gx[:, 2]-z0_gx)**2)
            ibu_gx = dbulge_gx <= R0_gx
            inlbulge_gx = snap_shot_bulge_gx[ibu_gx, :]
            
            if snap_N_gx[-2] != 0:
                dstar_gx = np.sqrt((snap_shot_star_gx[:, 0]-x0_gx)**2+(snap_shot_star_gx[:, 1]-y0_gx)**2 +(snap_shot_star_gx[:, 2]-z0_gx)**2)
                ids_gx = dstar_gx <= R0_gx
                inlstar_gx = snap_shot_star_gx[ids_gx, :]
                inlmass_star_gx = snap_mass_star_gx[ids_gx]
            if snap_N_gx[-1] != 0:
                dbndry_gx = np.sqrt((snap_shot_bndry_gx[:, 0]-x0_gx)**2+(snap_shot_bndry_gx[:, 1]-y0_gx)**2 +(snap_shot_bndry_gx[:, 2]-z0_gx)**2)
                idb_gx = dbndry_gx <= R0_gx
                inlbndry_gx = snap_shot_bndry_gx[idb_gx, :]
                inlmass_BH_gx = snap_mass_BH_gx[idb_gx]
            # try to find the most densitive area in the central area(which will set as BCG)
            R_range_gx = 100.0
            # set the central area as a cubic
            edge_x_gx = np.array([x0_gx-R_range_gx,x0_gx+R_range_gx])
            edge_y_gx = np.array([y0_gx-R_range_gx,y0_gx+R_range_gx])
            edge_z_gx = np.array([z0_gx-R_range_gx,z0_gx+R_range_gx])
            iv_s_gx = ((inlstar_gx[:,0]<=edge_x_gx[1])&(inlstar_gx[:,0]>=edge_x_gx[0]))&(
                       (inlstar_gx[:,1]<=edge_y_gx[1])&(inlstar_gx[:,1]>=edge_y_gx[0]))&((inlstar_gx[:,2]<=edge_z_gx[1])&(inlstar_gx[:,2]>=edge_z_gx[0]))
            inl_star_gx = inlstar_gx[iv_s_gx,:]
            # central postion and density of star
            num_bins_gx = np.ceil(R_range_gx*2/resolution)
            try :
                hist_star_gx, edge_star_gx = np.histogramdd(inl_star_gx, bins=(num_bins_gx, num_bins_gx, num_bins_gx))
            except ValueError:
                #print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z_gx)
                break
            bin_x_star_gx = np.array(edge_star_gx[0])
            bin_y_star_gx = np.array(edge_star_gx[1])
            bin_z_star_gx = np.array(edge_star_gx[2])
            inumber2 = hist_star_gx >= 10
            #next to find the center of BCG
            maxN_gx = len(inumber2)
            cen_po_star_gx = np.zeros((maxN_gx,3),dtype = np.float)
            hist_use_gx = hist_star_gx
            for q in range(maxN_gx):
                is_max_gx = np.unravel_index(np.argmax(hist_use_gx, axis=None), hist_use_gx.shape)
                cenxstar_gx = (bin_x_star_gx[is_max_gx[0]+1]+bin_x_star_gx[is_max_gx[0]])/2.0
                cenystar_gx = (bin_y_star_gx[is_max_gx[1]+1]+bin_y_star_gx[is_max_gx[1]])/2.0
                cenzstar_gx = (bin_z_star_gx[is_max_gx[2]+1]+bin_z_star_gx[is_max_gx[2]])/2.0
                cen_po_star_gx[q,:] = np.array([cenxstar_gx,cenystar_gx,cenzstar_gx])
                hist_use_gx[is_max_gx] = 0.0
            compare_d_gx= np.sqrt((cen_po_star_gx[:,0]-x0_gx)**2+
                                  (cen_po_star_gx[:,1]-y0_gx)**2+(cen_po_star_gx[:,2]-z0_gx)**2)
            ismin_gx = find.find1d(compare_d_gx,np.min(compare_d_gx))
            cen_x_star_gx = cen_po_star_gx[ismin_gx,0]
            cen_y_star_gx = cen_po_star_gx[ismin_gx,1]
            cen_z_star_gx = cen_po_star_gx[ismin_gx,2]
            # next,give radius and calculate the properties,in units kpc/h
            """
            according to Milky Way,assuming the BCG radiu is about 30kpc/h,but for finding the cut point
            give a 100kpc/h,to find the mean density profile
            """
            R_BCG_gx = size_BCG
            r_bcg_gx = np.logspace(1e-3,np.log10(R_BCG_gx),Nsize)
            n_r_gx = len(r_bcg_gx)
            inl_m_tot_gx = np.zeros(n_r_gx,dtype = np.float)
            inl_m_g_gx = np.zeros(n_r_gx,dtype = np.float)
            #above three mass instand:mass of a gas particle,star particle and DM particle aspectlly
            r_gas_gx = np.sqrt((inlgas_gx[:,0]-cen_x_star_gx)**2+(inlgas_gx[:,1]-cen_y_star_gx)**2+
                            (inlgas_gx[:,2]-cen_z_star_gx)**2)
            for q in range(n_r_gx):
                ia_gx = r_gas_gx <= r_bcg_gx[q]
                rgas_gx = r_gas_gx[ia_gx]
                inl_m_g_gx[q] = np.sum(inlmass_gas_gx[ia_gx]*10**10)
            m_g_gx[k] = inl_m_g_gx[-1]
            if snap_N_gx[-2] !=0:
                inl_m_s_gx = np.zeros(n_r_gx,dtype = np.float)
                r_star_gx = np.sqrt((inlstar_gx[:,0]-cen_x_star_gx)**2+(inlstar_gx[:,1]-cen_y_star_gx)**2+
                             (inlstar_gx[:,2]-cen_z_star_gx)**2)
                for q in range(n_r_gx):
                    ib_gx = r_star_gx <= r_bcg_gx[q]
                    rstar_gx = r_star_gx[ib_gx]
                    inl_m_s_gx[q] = np.sum(inlmass_star_gx[ib_gx]*10**10)
                m_s_gx[k] = inl_m_s_gx[-1]
            if snap_N_gx[-1] !=0:
                inl_BH_gx = {}
                inl_BH_m_gx = np.zeros(n_r_gx,dtype = np.float)
                r_BH_gx = np.sqrt((inlbndry_gx[:,0]-cen_x_star_gx)**2+(inlbndry_gx[:,1]-cen_y_star_gx)**2+
                             (inlbndry_gx[:,2]-cen_z_star_gx)**2)
                for q in range(n_r_gx):
                    ic_gx = r_BH_gx <=r_bcg_gx[q]
                    rbh_gx = r_BH_gx[ic_gx]
                    inl_BH_gx[q] = inlbndry_gx[ic_gx,:]
                    inl_BH_m_gx[q] = np.sum(inlmass_BH_gx[ic_gx]*10**10)
                m_bh_gx[k] = inl_BH_m_gx[-1]
            if snap_N_gx[-1] !=0 and snap_N_gx[-2] !=0:
                inl_m_tot_gx = inl_m_g_gx + inl_m_s_gx + inl_BH_m_gx
            if snap_N_gx[-1] ==0 and snap_N_gx[-2] !=0:
                inl_m_tot_gx = inl_m_g_gx + inl_m_s_gx
            if snap_N_gx[-1]==0 and snap_N_gx[-2] ==0:
                inl_m_tot_gx = inl_m_g_gx
            inl_m_b_gx = inl_m_g_gx + inl_m_s_gx
            inl_rho_b_gx = np.zeros(n_r_gx-1,dtype = np.float)
            inl_rho_g_gx = np.zeros(n_r_gx-1,dtype = np.float)
            inl_rho_s_gx = np.zeros(n_r_gx-1,dtype = np.float)
            for q in range(n_r_gx-1):
                dr_gx = r_bcg_gx[q+1] - r_bcg_gx[q]
                inl_rho_g_gx[q] = (inl_m_g_gx[q+1]-inl_m_g_gx[q])/(4*np.pi*dr_gx*r_bcg_gx[q]**2)
                if snap_N_gx[-2] !=0:
                    inl_rho_s_gx[q] = (inl_m_s_gx[q+1]-inl_m_s_gx[q])/(4*np.pi*dr_gx*r_bcg_gx[q]**2)
            rho_g_gx[k,:] = inl_rho_g_gx
            rho_s_gx[k,:] = inl_rho_s_gx
            inl_rho_b_gx = inl_rho_g_gx + inl_rho_s_gx
            mean_rho_g_gx = inl_m_g_gx/(4*np.pi*r_bcg_gx**3/3)
            m_rho_g_gx[k,:] = mean_rho_g_gx
            if snap_N_gx[-2] !=0:
                mean_rho_s_gx = inl_m_s_gx/(4*np.pi*r_bcg_gx**3/3)
                m_rho_s_gx[k,:] = mean_rho_s_gx
            mean_rho_b_gx = inl_m_b_gx/(4*np.pi*r_bcg_gx**3/3)      
    print('z_mu=',z_mu)
    print('z_gx=',z_gx)
    il = m_s != 0
    m_s = np.array(m_s[il])
    m_g = np.array(m_g[il])
    m_bh = np.array(m_bh[il])
    print('m_g=',m_g)
    print('m_s=',m_s)
    print('m_BH=',m_bh)
    m = np.array([m_s,m_g,m_bh])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_m_MU_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(m)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_m_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(np.int0(m.shape[0])):
            f['a'][t,:] = np.array(m[t,:])
    rho_s = np.array(rho_s[il,:])
    rho_g = np.array(rho_g[il,:])
    rho = np.array([rho_s,rho_g])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_rho_MU_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(rho)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_rho_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(np.int0(rho.shape[0])):
            f['a'][t,:] = np.array(rho[t,:])
    m_rho_s = np.array(m_rho_s[il,:])
    m_rho_g = np.array(m_rho_g[il,:])
    m_rho = np.array([m_rho_s,m_rho_g])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_mrho_MU_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(m_rho)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_mrho_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(np.int0(m_rho.shape[0])):
            f['a'][t,:] = np.array(m_rho[t,:])
    il_gx = m_s_gx != 0
    m_s_gx = np.array(m_s_gx[il_gx])
    m_g_gx = np.array(m_g_gx[il_gx])
    m_bh_gx = np.array(m_bh_gx[il_gx])
    print('m_g_gx=',m_g_gx)
    print('m_s_gx=',m_s_gx)
    print('m_BH_gx=',m_bh_gx)
    m_gx = np.array([m_s_gx,m_g_gx,m_bh_gx])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_m_GX_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(m_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_m_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(np.int0(m_gx.shape[0])):
            f['a'][t,:] = np.array(m_gx[t,:])
    rho_s_gx = np.array(rho_s_gx[il_gx,:])
    rho_g_gx = np.array(rho_g_gx[il_gx,:])
    rho_gx = np.array([rho_s_gx,rho_g_gx])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_rho_GX_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(rho_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_rho_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(np.int0(rho_gx.shape[0])):
            f['a'][t,:] = np.array(rho_gx[t,:])
    m_rho_s_gx = np.array(m_rho_s_gx[il_gx,:])
    m_rho_g_gx = np.array(m_rho_g_gx[il_gx,:])
    m_rho_gx = np.array([m_rho_s_gx,m_rho_g_gx])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_mrho_GX_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(m_rho_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_mrho_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(np.int0(m_rho_gx.shape[0])):
            f['a'][t,:] = np.array(m_rho_gx[t,:])
    z_gx = np.array(z_gx[il_gx])
    z_mu = np.array(z_mu[il])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_redshift_MU.h5','w') as f:
        f['a'] = np.array(z_mu)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_redshift_GX.h5','w') as f:
        f['a'] = np.array(z_gx)
    return
#BCG_evolution(resolution = True,g_size = True)
