"""
use for matter select
"""
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

resolution = 1.
_id_ = 0 # trace the first one halo(at z = 0)
BCG_size = 50 # in unit kpc/h
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)
R_in = 50.
R_out = 100.
def read_matter_MU():
    """
    read mass:
    m_BCG : (center distance) Cd <= 50 kpc/h (also use as mass of BCG)
    m_ICM : Cd > 100kpc/h but not in satellite
    md : dark matter mass
    m_sat : mass of satellite
    """
    z_mu = np.zeros(T_scale, dtype = np.float)
    m_d = np.zeros(T_scale, dtype = np.float)
    m_g = np.zeros(T_scale, dtype = np.float)
    m_s = np.zeros(T_scale, dtype = np.float)

    m_bcg = np.zeros(T_scale, dtype = np.float)
    m_icm = np.zeros(T_scale, dtype = np.float)
    m_sat = np.zeros(T_scale, dtype = np.float)

    m_icl = np.zeros(T_scale, dtype = np.float)
    mg_200 = np.zeros(T_scale, dtype = np.float)
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
            snap_shot_gas = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_mass_gas = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
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
            # total mass for each component
            try:
                x0 = xhalo[check_id]
                y0 = yhalo[check_id]
                z0 = zhalo[check_id]
                R0 = Rvir[check_id]
                m_g[k] = Mgas[check_id]
                m_s[k] = Mstar[check_id]
                m_d[k] = Mhalo[check_id]
            except IndexError:
                continue
            try:
                ih = sub_halo == halo_list[check_id]
                msat = Mstar[ih]
                sub_x = xhalo[ih]
                sub_y = yhalo[ih]
                sub_z = zhalo[ih]
                sub_r = Rvir[ih]

                dr_sub = np.sqrt((sub_x - x0)**2 
                    + (sub_y - y0)**2 + (sub_z - z0)**2)
                isub = dr_sub <= R0
                real_sub = msat[isub]
                real_sub_x = sub_x[isub]
                real_sub_y = sub_y[isub]
                real_sub_z = sub_z[isub]
                real_sub_r = sub_r[isub]

                mg_200[k] = np.sum(real_sub)
                cen_ip = np.where(real_sub == np.max(real_sub))[0]

                cen_mass = real_sub[cen_ip[0]]
                m_sat[k] = np.sum(real_sub) - cen_mass

                cen_x = real_sub_x[cen_ip[0]]
                cen_y = real_sub_y[cen_ip[0]]
                cen_z = real_sub_z[cen_ip[0]]

                # for star 
                drs = np.sqrt((snap_shot_star[:, 0] - cen_x)**2 +
                              (snap_shot_star[:, 1] - cen_y)**2 +(snap_shot_star[:, 2] - cen_z)**2)

                isr1 = drs <= R_in
                inlstar1 = snap_shot_star[isr1, :]
                inlmass_star1 = snap_mass_star[isr1]
                m_bcg[k] = np.sum(snap_mass_star[isr1])*10**10
                
                isr3 = (drs > R_out) & (drs <= R0)
                inlstar3 = snap_shot_star[isr3, :]
                inlmass_star3 = snap_mass_star[isr3]
                
               	m_tt = 0
               	pos = inlstar3*1
               	out_star = inlmass_star3*1
                for q in range(len(sub_r)):
                    dr = np.sqrt((pos[:,0] - sub_x[q])**2 + 
                        (pos[:,1] - sub_y[q])**2 + (pos[:,2] - sub_z[q])**2)

                    iic = dr >= sub_r[q]
                    pos = pos[iic, :]
                    out_star = out_star[iic]
                m_icm[k] = np.sum(out_star)*10**10

                m_icl[k] = m_s[k] - m_bcg[k] - m_sat[k]
            except ValueError:
                m_g[k] = Mgas[check_id]
                m_s[k] = Mstar[check_id]
                m_d[k] = Mhalo[check_id]
                m_sat[k] = 0
                m_icm[k] = 0
                mg_200[k] = 0
                m_bcg[k] = Mstar[check_id]

                m_icl[k] = m_s[k] - m_bcg[k] - m_sat[k]
                continue

    il = m_s != 0
    ms = np.array(m_s[il])
    mg = np.array(m_g[il])
    md = np.array(m_d[il])
    zmu = np.array(z_mu[il])

    mbcg = np.array(m_bcg[il])
    micm = np.array(m_icm[il])
    msat = np.array(m_sat[il])

    micl = np.array(m_icl[il])
    mg200 = np.array(mg_200[il])
    data = np.array([zmu, ms, mg, md, mbcg, micm, msat, micl, mg200])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/MUSIC_matter_select_bason_halo_Cen%.1f.h5'%BCG_size, 'w') as f:
        f['a'] = np.array(data)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/MUSIC_matter_select_bason_halo_Cen%.1f.h5'%BCG_size) as f:
        for q in range(len(data)):
            f['a'][q, :] = data[q, :]

    return

def read_matter_GX():

    z_gx = np.zeros(T_scale, dtype = np.float)
    m_d = np.zeros(T_scale, dtype = np.float)
    m_g = np.zeros(T_scale, dtype = np.float)
    m_s = np.zeros(T_scale, dtype = np.float)

    m_bcg = np.zeros(T_scale, dtype = np.float)
    m_icm = np.zeros(T_scale, dtype = np.float)
    m_sat = np.zeros(T_scale, dtype = np.float)

    m_icl = np.zeros(T_scale, dtype = np.float)
    mg_200 = np.zeros(T_scale, dtype = np.float)
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
        print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
    	      % (snap_N[0], snap_N[1], snap_N[2], snap_N[3], snap_N[4], snap_N[5]))
        if snap_z <= id_z:
            z_gx[k] = snap_z
            snap_shot_gas = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_mass_gas = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
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
            goalz = find.find1d(com_z,goal_z)
            goal_halo = tree_line[_id_,goalz]
            halo_list = np.array(main_halo['col1'])
            sub_halo = np.array(main_halo['col2'])

            Rvir = np.array(main_halo['col12'])
            xhalo = np.array(main_halo['col6'])
            yhalo = np.array(main_halo['col7'])
            zhalo = np.array(main_halo['col8'])
            Mhalo = np.array(main_halo['col4'])
            Mgas = np.array(main_halo['col45'])
            Mstar= np.array(main_halo['col65'])

            check_id = find.find1d(halo_list, goal_halo)
            # total mass for each component
            try:
                x0 = xhalo[check_id]
                y0 = yhalo[check_id]
                z0 = zhalo[check_id]
                R0 = Rvir[check_id]
                m_g[k] = Mgas[check_id]
                m_s[k] = Mstar[check_id]
                m_d[k] = Mhalo[check_id]
            except IndexError:
                continue
            try:
                ih = sub_halo == halo_list[check_id]
                msat = Mstar[ih]
                sub_x = xhalo[ih]
                sub_y = yhalo[ih]
                sub_z = zhalo[ih]
                sub_r = Rvir[ih]

                dr_sub = np.sqrt((sub_x - x0)**2 
                    + (sub_y - y0)**2 + (sub_z - z0)**2)
                isub = dr_sub <= R0
                real_sub = msat[isub]
                real_sub_x = sub_x[isub]
                real_sub_y = sub_y[isub]
                real_sub_z = sub_z[isub]
                real_sub_r = sub_r[isub]

                mg_200[k] = np.sum(real_sub)
                cen_ip = np.where(real_sub == np.max(real_sub))[0]
                
                cen_mass = real_sub[cen_ip[0]]
                m_sat[k] = np.sum(real_sub) - cen_mass

                cen_x = real_sub_x[cen_ip[0]]
                cen_y = real_sub_y[cen_ip[0]]
                cen_z = real_sub_z[cen_ip[0]]

                # for star 
                drs = np.sqrt((snap_shot_star[:, 0] - cen_x)**2 +
                              (snap_shot_star[:, 1] - cen_y)**2 +(snap_shot_star[:, 2] - cen_z)**2)

                isr1 = drs <= R_in
                inlstar1 = snap_shot_star[isr1, :]
                inlmass_star1 = snap_mass_star[isr1]
                m_bcg[k] = np.sum(snap_mass_star[isr1])*10**10

                isr3 = (drs > R_out) & (drs <= R0)
                inlstar3 = snap_shot_star[isr3, :]
                inlmass_star3 = snap_mass_star[isr3]

                m_tt = 0
                pos = inlstar3*1
                out_star = inlmass_star3*1
                for q in range(len(sub_r)):
                    dr = np.sqrt((pos[:,0] - sub_x[q])**2 + 
                        (pos[:,1] - sub_y[q])**2 + (pos[:,2] - sub_z[q])**2)

                    iic = dr >= sub_r[q]
                    pos = pos[iic, :]
                    out_star = out_star[iic]
                m_icm[k] = np.sum(out_star)*10**10

                m_icl[k] = m_s[k] - m_bcg[k] - m_sat[k]
            except ValueError:
                m_g[k] = Mgas[check_id]
                m_s[k] = Mstar[check_id]
                m_d[k] = Mhalo[check_id]
                m_sat[k] = 0
                m_icm[k] = 0
                mg_200[k] = 0
                m_bcg[k] = Mstar[check_id]

                m_icl[k] = m_s[k] - m_bcg[k] - m_sat[k]
                continue

    il = m_s != 0
    ms = np.array(m_s[il])
    mg = np.array(m_g[il])
    md = np.array(m_d[il])
    zgx = np.array(z_gx[il])

    mbcg = np.array(m_bcg[il])
    micm = np.array(m_icm[il])
    msat = np.array(m_sat[il])

    micl = np.array(m_icl[il])
    mg200 = np.array(mg_200[il])
    data = np.array([zgx, ms, mg, md, mbcg, micm, msat, micl, mg200])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/GadgetX_matter_select_bason_halo_Cen%.1f.h5'%BCG_size, 'w') as f:
        f['a'] = np.array(data)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/GadgetX_matter_select_bason_halo_Cen%.1f.h5'%BCG_size) as f:
        for q in range(len(data)):
            f['a'][q, :] = data[q, :]


def fig_out():
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/GadgetX_matter_select_bason_halo_Cen%.1f.h5'%BCG_size) as f:
        data_x = np.array(f['a'])
    z_gx = data_x[0,:]
    ms_gx = data_x[1,:]
    mg_gx = data_x[2,:]
    md_gx = data_x[3,:]
    mbcg_gx = data_x[4,:]
    micm_gx = data_x[5,:]
    msat_gx = data_x[6,:]

    ix = msat_gx != 0
    zgx = z_gx[ix]
    msgx = ms_gx[ix]
    mggx = mg_gx[ix]
    mdgx = md_gx[ix]
    mbcggx = mbcg_gx[ix]
    micmgx = micm_gx[ix]
    msatgx = msat_gx[ix]

    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/MUSIC_matter_select_bason_halo_Cen%.1f.h5'%BCG_size) as f:
        data_m = np.array(f['a'])
    z_mu = data_m[0,:]
    ms_mu = data_m[1,:]
    mg_mu = data_m[2,:]
    md_mu = data_m[3,:]
    mbcg_mu = data_m[4,:]
    micm_mu = data_m[5,:]
    msat_mu = data_m[6,:]

    iy = msat_mu != 0
    zmu = z_mu[iy]
    msmu = ms_mu[iy]
    mgmu = mg_mu[iy]
    mdmu = md_mu[iy]
    mbcgmu = mbcg_mu[iy]
    micmmu = micm_mu[iy]
    msatmu = msat_mu[iy]

    plt.figure()
    plt.plot(zgx, mdgx, ls = '-', c = 'k', label = r'$M_{DM}$', )
    plt.plot(zgx, mbcggx, ls = '-', c = 'r', label = r'$M_{BCG}$',)
    plt.plot(zgx, micmgx, ls = '-', c = 'g', label = r'$M_{ICM}$',)
    plt.plot(zgx, msatgx, ls = '-', c = 'b', label = r'$M_{satellite}$')
    plt.tick_params(axis = 'both', which = 'both', direction = 'in')
    plt.title('Mass evolution in GadgetX')
    plt.xlabel('z')
    plt.ylabel(r'$M[M_\odot /h]$')
    plt.yscale('log')
    plt.legend(loc = 1)
    plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/mass_evolution_with_ICM_GX.png', dpi = 300)
    plt.close()

    plt.figure()
    plt.plot(zmu, mdmu, ls = '-', c = 'k', label = r'$M_{DM}$',)
    plt.plot(zmu, mbcgmu, ls = '-', c = 'r', label = r'$M_{BCG}$',)
    plt.plot(zmu, micmmu, ls = '-', c = 'g', label = r'$M_{ICM}$',)
    plt.plot(zmu, msatmu, ls = '-', c = 'b', label = r'$M_{satellite}$',)
    plt.tick_params(axis = 'both', which = 'both', direction = 'in')
    plt.legend(loc = 1)
    plt.title('Mass evolution in MUSIC')
    plt.xlabel('z')
    plt.ylabel(r'$M[M_\odot /h]$')
    plt.yscale('log')
    plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/mass_evolution_with_ICM_MU.png', dpi = 300)
    plt.close()

    return
def main():
    read_matter_MU()
    read_matter_GX()
    fig_out()
if __name__ == "__main__":
    main()