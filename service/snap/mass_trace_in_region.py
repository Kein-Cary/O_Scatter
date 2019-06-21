import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np
import astropy.io.ascii as asc
import h5py

import find
import changds
import pygadgetreader as pygdr
#import handy.scatter as hsc
"""
ip : the goal halo number of halo at z==0
resolution : the resolution of graph
viewp : the view direction,must be integer
# view control the projection of distribution：1--x-y panel,2--x-z panel,3--y-z panel
"""
h = 0.7
resolution = 1. # in unit kpc/h
_id_ = 0 # trace the first one halo(at z = 0)
BCG_size = 50 # in unit kpc/h
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)
R_in = 50
R_out = 100
def read_data_MUSIC():

    m_s = np.zeros(T_scale, dtype = np.float)
    m_d = np.zeros(T_scale, dtype = np.float)

    m_sat = np.zeros(T_scale, dtype = np.float)
    m_bcg = np.zeros(T_scale, dtype = np.float)
    m_ICL = np.zeros(T_scale, dtype = np.float)
    m_ICM = np.zeros(T_scale, dtype = np.float)

    m_s_iner = np.zeros(T_scale, dtype = np.float)
    m_s_media = np.zeros(T_scale, dtype = np.float)
    m_s_out = np.zeros(T_scale, dtype = np.float)

    z_mu = np.zeros(T_scale,dtype = np.float)
    # read the meeger tree
    with h5py.File('/home/cxkttwl/Scatter/MUSIC/Redshift.h5') as f:
        y0 = f['a']
        com_z = np.array(y0)
    with h5py.File('/home/cxkttwl/Scatter/MUSIC/main_tree.h5') as f:
        y1 = f['a']
        tree_line = np.array(y1)
    L = tree_line.shape[0]
    iv = np.zeros(L,dtype = np.int0)
    for l in range(L):
        u = find.find1d(tree_line[l,:],0)
        iv[l] = u
    id_z = com_z[iv[_id_]]

    for k in range(T_scale):
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
        if snap_z <= id_z:
            z_mu[k] = snap_z
            try:
                snap_shot_star = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
                snap_mass_star = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
            except SystemExit:
                print('no star particles now')
                break
            main_halo = asc.read(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z,
                                 converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z = changds.inv_chidas('%.3f'%snap_z)
            goalz = find.find1d(com_z,goal_z)
            goal_halo = tree_line[_id_,goalz]
            halo_list = np.array(main_halo['col1'])
            sub_halo = np.array(main_halo['col2'])
            check_id = find.find1d(halo_list, goal_halo)
            Rvir = np.array(main_halo['col12'])
            xhalo = np.array(main_halo['col6'])
            yhalo = np.array(main_halo['col7'])
            zhalo = np.array(main_halo['col8'])
            Mhalo = np.array(main_halo['col4'])
            Mgas = np.array(main_halo['col45'])
            Mstar= np.array(main_halo['col65'])
            # select the total star, gas, dark matter mass
            try:
                x0 = xhalo[check_id]
                y0 = yhalo[check_id]
                z0 = zhalo[check_id]
                R0 = Rvir[check_id]
                m_d[k] = Mhalo[check_id]

                if snap_N[-2] != 0:
                    dstar = np.sqrt((snap_shot_star[:, 0] - x0)**2 + 
                        (snap_shot_star[:, 1] - y0)**2 + (snap_shot_star[:, 2] - z0)**2)
                    ids = dstar <= R0
                    inlstar = snap_shot_star[ids, :]
                    inlmass_star = snap_mass_star[ids]
                m_s[k] = np.sum(inlmass_star) * 10**10

                R_range = 150.0
                edge_x = np.array([x0-R_range, x0+R_range])
                edge_y = np.array([y0-R_range, y0+R_range])
                edge_z = np.array([z0-R_range, z0+R_range])
                iv_s = ((inlstar[:,0] <= edge_x[1])&(inlstar[:,0] >= edge_x[0]))&(
                        (inlstar[:,1] <= edge_y[1])&(inlstar[:,1] >= edge_y[0]))&(
                        (inlstar[:,2] <= edge_z[1])&(inlstar[:,2] >= edge_z[0]))
                inl_star = inlstar[iv_s,:]
                num_bins = np.ceil(R_range*2/resolution)
                try:
                    hist_star, edge_star = np.histogramdd(inl_star, bins=(num_bins, num_bins, num_bins))                
                    bin_x_star = np.array(edge_star[0])
                    bin_y_star = np.array(edge_star[1])
                    bin_z_star = np.array(edge_star[2])

                    inumber1 = hist_star >= 10
                    maxN = len(inumber1)
                    cen_po_star = np.zeros((maxN,3),dtype = np.float)
                    hist_use = hist_star
                    for p in range(maxN):
                        is_max = np.unravel_index(np.argmax(hist_use, axis=None), hist_use.shape)
                        cenxstar = (bin_x_star[is_max[0]+1] + bin_x_star[is_max[0]])/2.0
                        cenystar = (bin_y_star[is_max[1]+1] + bin_y_star[is_max[1]])/2.0
                        cenzstar = (bin_z_star[is_max[2]+1] + bin_z_star[is_max[2]])/2.0
                        cen_po_star[p,:] = np.array([cenxstar,cenystar,cenzstar])
                        hist_use[is_max] = 0.0
                    compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+
                                       (cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
                    ismin = find.find1d(compare_d,np.min(compare_d))
                    cen_x_star = cen_po_star[ismin,0]
                    cen_y_star = cen_po_star[ismin,1]
                    cen_z_star = cen_po_star[ismin,2]
                    # calculate the bcg mass
                    ds_star = np.sqrt((inlstar[:,0] - cen_x_star)**2 + 
                        (inlstar[:,1] - cen_y_star)**2 + (inlstar[:,2] - cen_z_star)**2)
                    icen = ds_star <= BCG_size
                    gx = np.where(icen == True)[0]
                    m_bcg[k] = np.sum(inlmass_star[gx]) * 10**10
                    
                    # select the satellite and ICL
                    try:
                        ih = sub_halo == halo_list[check_id]
                        msat = Mstar[ih]
                        sub_x = xhalo[ih]
                        sub_y = yhalo[ih]
                        sub_z = zhalo[ih]
                        sub_r = Rvir[ih]

                        dr_sub0 = np.sqrt((sub_x - x0)**2 
                            + (sub_y - y0)**2 + (sub_z - z0)**2)
                        isub = dr_sub0 <= R0
                        real_sub = msat[isub]
                        real_sub_x = sub_x[isub]
                        real_sub_y = sub_y[isub]
                        real_sub_z = sub_z[isub]
                        real_sub_r = sub_r[isub]

                        dr_sub1 = np.sqrt((real_sub_x - cen_x_star)**2 + 
                            (real_sub_y - cen_y_star)**2 + (real_sub_z - cen_z_star)**2)
                        ucen = np.where(dr_sub1 == np.min(dr_sub1))[0]
                        cen_mass = real_sub[ucen[0]]
                        m_sat[k] = np.sum(real_sub) - cen_mass
                        m_ICL[k] = m_s[k] - m_sat[k] - m_bcg[k]
                    except ValueError:
                        print('there is noly one halo!')
                        print('Now redshift is %.3f'%snap_z)
                        m_sat[k] = 0
                        m_ICL[k] = m_s[k] - m_sat[k] - m_bcg[k]

                    # for star 
                    drs = np.sqrt((inlstar[:, 0] - cen_x_star)**2 +
                                  (inlstar[:, 1] - cen_y_star)**2 +(inlstar[:, 2] - cen_z_star)**2)
                    isr1 = drs <= R_in
                    inlstar1 = inlstar[isr1, :]
                    inlmass_star1 = inlmass_star[isr1]
                    m_s_iner[k] = np.sum(inlmass_star[isr1])*10**10
                    
                    isr2 = (drs > R_in) & (drs <= R_out)
                    inlstar2 = inlstar[isr2, :]
                    inlmass_star2 = inlmass_star[isr2]
                    m_s_media[k] = np.sum(inlmass_star[isr2])*10**10
                    
                    isr3 = (drs > R_out) & (drs <= R0)
                    inlstar3 = inlstar[isr3, :]
                    inlmass_star3 = inlmass_star[isr3]
                    m_s_out[k] = np.sum(inlmass_star[isr3])*10**10

                    pos = inlstar3*1
                    out_star = inlmass_star3*1
                    for q in range(len(sub_r)):
                        dr = np.sqrt((pos[:,0] - sub_x[q])**2 + 
                            (pos[:,1] - sub_y[q])**2 + (pos[:,2] - sub_z[q])**2)

                        iic = dr >= sub_r[q]
                        pos = pos[iic, :]
                        out_star = out_star[iic]
                    m_ICM[k] = np.sum(out_star)*10**10
                except ValueError:
                    print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
                    break
            except IndexError:
                print('There is no halo to trace!')
                break
    ### save data
    il = m_s != 0
    zmu = np.array(z_mu[il])
    ms = np.array(m_s[il])
    md = np.array(m_d[il])

    ms_iner = np.array(m_s_iner[il])
    ms_media = np.array(m_s_media[il])
    ms_out = np.array(m_s_out[il])
    
    ms_sat = np.array(m_sat[il])
    ms_ICL = np.array(m_ICL[il])
    ms_bcg = np.array(m_bcg[il])
    ms_ICM = np.array(m_ICM[il])
    region_data = np.array([zmu, ms, md, ms_iner, ms_media, ms_out, ms_sat, ms_ICL, ms_bcg, ms_ICM])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/MUSIC_trace_in_region_BCG%.1f_Reso%.3f.h5'%(BCG_size, resolution), 'w') as f:
        f['a'] = np.array(region_data)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/MUSIC_trace_in_region_BCG%.1f_Reso%.3f.h5'%(BCG_size, resolution)) as f:
        for q in range(len(region_data)):
            f['a'][q, :] = region_data[q, :]

    return

def read_data_GadgetX():

    m_s = np.zeros(T_scale, dtype = np.float)
    m_d = np.zeros(T_scale, dtype = np.float)

    m_sat = np.zeros(T_scale, dtype = np.float)
    m_bcg = np.zeros(T_scale, dtype = np.float)
    m_ICL = np.zeros(T_scale, dtype = np.float)
    m_ICM = np.zeros(T_scale, dtype = np.float)

    m_s_iner = np.zeros(T_scale, dtype = np.float)
    m_s_media = np.zeros(T_scale, dtype = np.float)
    m_s_out = np.zeros(T_scale, dtype = np.float)

    z_gx = np.zeros(T_scale,dtype = np.float)
    with h5py.File('/home/cxkttwl/Scatter/GadgetX/Redshift_GX.h5') as f:
        y2 = f['a']
        com_z = np.array(y2)
    with h5py.File('/home/cxkttwl/Scatter/GadgetX/main_tree_GX.h5') as f:
        y3 = f['a']
        tree_line = np.array(y3)
    L = tree_line.shape[0]
    iv = np.zeros(L,dtype = np.int0)
    for ll in range(L):
        u = find.find1d(tree_line[ll,:],0)
        iv[ll] = u
    id_z = com_z[iv[_id_]]

    for k in range(T_scale):
        dd = alpha[-k]
        ss = str(dd)
        if len(ss) == 1:
            s_u = str('00%.0f' % dd)
        elif len(ss) == 2:
            s_u = str('0%.0f' % dd)
        else:
            s_u = str(dd)
        No_snap = s_u
        snap_z = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
        snap_z = np.abs(snap_z)
        print('Now redshift is %.3f' % snap_z)
        snap_N = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
        if snap_z <= id_z:
            z_gx[k] = snap_z
            try:
                snap_shot_star = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
                snap_mass_star = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
            except SystemExit:
                print('no star particles now')
                break
            main_halo = asc.read(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z,
                                converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z = changds.inv_chidas('%.3f'%snap_z)
            goalz = find.find1d(com_z, goal_z)
            goal_halo = tree_line[_id_, goalz]
            halo_list = np.array(main_halo['col1'])
            sub_halo = np.array(main_halo['col2'])
            check_id = find.find1d(halo_list, goal_halo)
            Rvir = np.array(main_halo['col12'])
            xhalo = np.array(main_halo['col6'])
            yhalo = np.array(main_halo['col7'])
            zhalo = np.array(main_halo['col8'])
            Mhalo = np.array(main_halo['col4'])
            Mgas = np.array(main_halo['col45'])
            Mstar= np.array(main_halo['col65'])
            # select the total star, gas, dark matter mass
            try:
                x0 = xhalo[check_id]
                y0 = yhalo[check_id]
                z0 = zhalo[check_id]
                R0 = Rvir[check_id]
                m_d[k] = Mhalo[check_id]

                if snap_N[-2] != 0:
                    dstar = np.sqrt((snap_shot_star[:, 0] - x0)**2 + 
                        (snap_shot_star[:, 1] - y0)**2 + (snap_shot_star[:, 2] - z0)**2)
                    ids = dstar <= R0
                    inlstar = snap_shot_star[ids, :]
                    inlmass_star = snap_mass_star[ids]
                m_s[k] = np.sum(inlmass_star) * 10**10

                R_range = 150.0
                edge_x = np.array([x0-R_range, x0+R_range])
                edge_y = np.array([y0-R_range, y0+R_range])
                edge_z = np.array([z0-R_range, z0+R_range])
                iv_s = ((inlstar[:,0] <= edge_x[1])&(inlstar[:,0] >= edge_x[0]))&(
                        (inlstar[:,1] <= edge_y[1])&(inlstar[:,1] >= edge_y[0]))&(
                        (inlstar[:,2] <= edge_z[1])&(inlstar[:,2] >= edge_z[0]))
                inl_star = inlstar[iv_s,:]
                num_bins = np.ceil(R_range*2/resolution)
                try:
                    hist_star, edge_star = np.histogramdd(inl_star, bins=(num_bins, num_bins, num_bins))
                    bin_x_star = np.array(edge_star[0])
                    bin_y_star = np.array(edge_star[1])
                    bin_z_star = np.array(edge_star[2])

                    inumber1 = hist_star >= 10
                    maxN = len(inumber1)
                    cen_po_star = np.zeros((maxN,3),dtype = np.float)
                    hist_use = hist_star
                    for p in range(maxN):
                        is_max = np.unravel_index(np.argmax(hist_use, axis=None), hist_use.shape)
                        cenxstar = (bin_x_star[is_max[0]+1] + bin_x_star[is_max[0]])/2.0
                        cenystar = (bin_y_star[is_max[1]+1] + bin_y_star[is_max[1]])/2.0
                        cenzstar = (bin_z_star[is_max[2]+1] + bin_z_star[is_max[2]])/2.0
                        cen_po_star[p,:] = np.array([cenxstar,cenystar,cenzstar])
                        hist_use[is_max] = 0.0
                    compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+
                                       (cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
                    ismin = find.find1d(compare_d,np.min(compare_d))
                    cen_x_star = cen_po_star[ismin,0]
                    cen_y_star = cen_po_star[ismin,1]
                    cen_z_star = cen_po_star[ismin,2]
                    # calculate the bcg mass
                    ds_star = np.sqrt((inlstar[:,0] - cen_x_star)**2 + 
                        (inlstar[:,1] - cen_y_star)**2 + (inlstar[:,2] - cen_z_star)**2)
                    icen = ds_star <= BCG_size
                    gx = np.where(icen == True)[0]
                    m_bcg[k] = np.sum(inlmass_star[gx]) * 10**10

                    # select the satellite and ICL
                    try:
                        ih = sub_halo == halo_list[check_id]
                        msat = Mstar[ih]
                        sub_x = xhalo[ih]
                        sub_y = yhalo[ih]
                        sub_z = zhalo[ih]
                        sub_r = Rvir[ih]

                        dr_sub0 = np.sqrt((sub_x - x0)**2 
                            + (sub_y - y0)**2 + (sub_z - z0)**2)
                        isub = dr_sub0 <= R0
                        real_sub = msat[isub]
                        real_sub_x = sub_x[isub]
                        real_sub_y = sub_y[isub]
                        real_sub_z = sub_z[isub]
                        real_sub_r = sub_r[isub]

                        dr_sub1 = np.sqrt((real_sub_x - cen_x_star)**2 + 
                            (real_sub_y - cen_y_star)**2 + (real_sub_z - cen_z_star)**2)
                        ucen = np.where(dr_sub1 == np.min(dr_sub1))[0]
                        cen_mass = real_sub[ucen[0]]
                        m_sat[k] = np.sum(real_sub) - cen_mass
                        m_ICL[k] = m_s[k] - m_sat[k] - m_bcg[k]
                    except ValueError:
                        print('there is noly one halo!')
                        print('Now redshift is %.3f'%snap_z)
                        m_sat[k] = 0
                        m_ICL[k] = m_s[k] - m_sat[k] - m_bcg[k]

                    # for star 
                    drs = np.sqrt((inlstar[:, 0] - cen_x_star)**2 +
                                  (inlstar[:, 1] - cen_y_star)**2 +(inlstar[:, 2] - cen_z_star)**2)
                    isr1 = drs <= R_in
                    inlstar1 = inlstar[isr1, :]
                    inlmass_star1 = inlmass_star[isr1]
                    m_s_iner[k] = np.sum(inlmass_star[isr1])*10**10
                    
                    isr2 = (drs > R_in) & (drs <= R_out)
                    inlstar2 = inlstar[isr2, :]
                    inlmass_star2 = inlmass_star[isr2]
                    m_s_media[k] = np.sum(inlmass_star[isr2])*10**10
                    
                    isr3 = (drs > R_out) & (drs <= R0)
                    inlstar3 = inlstar[isr3, :]
                    inlmass_star3 = inlmass_star[isr3]
                    m_s_out[k] = np.sum(inlmass_star[isr3])*10**10

                    pos = inlstar3*1
                    out_star = inlmass_star3*1
                    for q in range(len(sub_r)):
                        dr = np.sqrt((pos[:,0] - sub_x[q])**2 + 
                            (pos[:,1] - sub_y[q])**2 + (pos[:,2] - sub_z[q])**2)

                        iic = dr >= sub_r[q]
                        pos = pos[iic, :]
                        out_star = out_star[iic]
                    m_ICM[k] = np.sum(out_star)*10**10
                except ValueError:
                    print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
                    break
            except IndexError:
                print('There is no halo to trace!')
                break
    ### save data
    il = m_s != 0
    zgx = np.array(z_gx[il])
    ms = np.array(m_s[il])
    md = np.array(m_d[il])

    ms_iner = np.array(m_s_iner[il])
    ms_media = np.array(m_s_media[il])
    ms_out = np.array(m_s_out[il])
    
    ms_sat = np.array(m_sat[il])
    ms_ICL = np.array(m_ICL[il])
    ms_bcg = np.array(m_bcg[il])
    ms_ICM = np.array(m_ICM[il])
    region_data = np.array([zgx, ms, md, 
                ms_iner, ms_media, ms_out, 
                ms_sat, ms_ICL, ms_bcg, ms_ICM])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/GadgetX_trace_in_region_BCG%.1f_Reso%.3f.h5'%(BCG_size, resolution), 'w') as f:
        f['a'] = np.array(region_data)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/GadgetX_trace_in_region_BCG%.1f_Reso%.3f.h5'%(BCG_size, resolution)) as f:
        for q in range(len(region_data)):
            f['a'][q, :] = region_data[q, :]

    return

def fig_result():
    
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/GadgetX_trace_in_region_BCG%.1f_Reso%.3f.h5'%(BCG_size, resolution)) as f:
        GadgetX = np.array(f['a'])
    z_gx = GadgetX[0,:]
    ms_gx = GadgetX[1,:]
    md_gx = GadgetX[2,:]
    
    ms_iner_gx = GadgetX[3,:]
    ms_media_gx = GadgetX[4,:]
    ms_out_gx = GadgetX[5,:]
    
    m_sat_gx = GadgetX[6,:]
    m_ICL_gx = GadgetX[7,:]
    m_bcg_gx = GadgetX[8,:]
    m_ICM_gx = GadgetX[9,:]

    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/MUSIC_trace_in_region_BCG%.1f_Reso%.3f.h5'%(BCG_size, resolution)) as f:
        Music = np.array(f['a'])
    z_mu = Music[0,:]
    ms_mu = Music[1,:]
    md_mu = Music[2,:]
    
    ms_iner_mu = Music[3,:]
    ms_media_mu = Music[4,:]
    ms_out_mu = Music[5,:]
    
    m_sat_mu = Music[6,:]
    m_ICL_mu = Music[7,:]
    m_bcg_mu = Music[8,:]
    m_ICM_mu = Music[9,:]

    plt.figure()
    plt.plot(z_mu[md_mu != 0], md_mu[md_mu != 0], 'k-', label = r'$M_{DM}$')
    plt.plot(z_mu[m_sat_mu != 0], m_sat_mu[m_sat_mu != 0], 'g-', label = r'$M_{satellite}$')
    plt.plot(z_mu[m_bcg_mu != 0], m_bcg_mu[m_bcg_mu != 0], 'r-', label = r'$M_{BCG}$')
    plt.plot(z_mu[m_ICM_mu != 0], m_ICM_mu[m_ICM_mu != 0], 'b-', label = r'$M_{ICM}$')
    plt.legend(loc = 1)
    plt.tick_params(axis = 'both', which = 'both', direction = 'in')
    plt.xlabel('z')
    plt.ylabel(r'$Mass[M_\odot /h]$')
    plt.yscale('log')
    plt.title('Mass as function of z in MUSIC')
    plt.xlim(0, 5)
    plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/MUSIC_mass_evolution_bas_on_star.png', dpi = 600)
    plt.close()

    plt.figure()
    plt.plot(z_gx[md_gx != 0], md_gx[md_gx != 0], 'k-', label = r'$M_{DM}$')
    plt.plot(z_gx[m_sat_gx != 0], m_sat_gx[m_sat_gx != 0], 'g-', label = r'$M_{satellite}$')
    plt.plot(z_gx[m_bcg_gx != 0], m_bcg_gx[m_bcg_gx != 0], 'r-', label = r'$M_{BCG}$')
    plt.plot(z_gx[m_ICM_gx != 0], m_ICM_gx[m_ICM_gx != 0], 'b-', label = r'$M_{ICM}$')
    plt.legend(loc = 1)
    plt.tick_params(axis = 'both', which = 'both', direction = 'in')
    plt.xlabel('z')
    plt.ylabel(r'$Mass[M_\odot /h]$')
    plt.yscale('log')
    plt.title('Mass as function of z in GadgetX')
    plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/GadgetX_mass_evolution_bas_on_star.png', dpi = 600)
    plt.close()

    plt.figure()
    plt.plot(z_gx[md_gx != 0], md_gx[md_gx != 0], ls = '-', c = 'k', label = r'$M_{DM}$')
    plt.plot(z_gx[ms_iner_gx != 0], ms_iner_gx[ms_iner_gx != 0], ls = '-', c = 'r', label = r'$[M_\ast]_{R \leq 50kpc/h}$')
    plt.plot(z_gx[ms_media_gx != 0], ms_media_gx[ms_media_gx != 0], ls = '-', c = 'b', label = r'$[M_\ast]_{50kpc/h < R \leq 100kpc/h}$')
    plt.plot(z_gx[ms_out_gx != 0], ms_out_gx[ms_out_gx != 0], ls = '-', c = 'g', label = r'$[M_\ast]_{100kpc/h \leq R}$')
    plt.tick_params(axis = 'both', which = 'both', direction = 'in')
    plt.legend(loc = 1)
    plt.xlabel('z')
    plt.ylabel(r'$Mass[M_\odot/h]$')
    plt.yscale('log')
    plt.title('radiu mass distribution in GadgetX')
    plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/mass_contribut_GadgetX_bas_on_star.png', dpi = 600)
    plt.close()  

    plt.figure()
    plt.plot(z_mu[md_mu != 0], md_mu[md_mu != 0], 'k-', label = r'$M_{DM}$')
    plt.plot(z_mu[ms_iner_mu != 0], ms_iner_mu[ms_iner_mu != 0], 'r-', label = r'$[M_\ast]_{R \leq 50kpc/h}$')
    plt.plot(z_mu[ms_media_mu != 0], ms_media_mu[ms_media_mu != 0], 'g-', label = r'$[M_\ast]_{50kpc/h < R \leq 100kpc/h}$')
    plt.plot(z_mu[ms_out_mu != 0], ms_out_mu[ms_out_mu != 0], 'b-', label= r'$[M_\ast]_{100kpc/h \leq R}$')
    plt.legend(loc = 1)
    plt.tick_params(axis = 'both', which = 'both', direction = 'in')
    plt.xlabel('z')
    plt.ylabel(r'$Mass[M_\odot /h]$')
    plt.yscale('log')
    plt.title('radiu mass distribution in MUSIC')
    plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/mass_contribut_MUSIC_bas_on_star.png', dpi = 600)
    plt.close()
    raise
    return

def main():
    #read_data_MUSIC()
    #read_data_GadgetX()
    fig_result()
    
if __name__ == "__main__":
    main()