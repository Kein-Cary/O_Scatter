# this file try to find the BCG in a given halo
import matplotlib as mpl
mpl.use('Agg')
#from pygadgetreader import *
import pygadgetreader as pygdr
import numpy as np
import h5py
import astropy.io.ascii as asc
import matplotlib.pyplot as plt
#from handy import scatter as hsc
import find
import changds 
import round_number as rnr
def find_BCG(ip,resolution,g_size):
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
    #save the center of BCG
    Read_cen_gla_MU = np.zeros((len(alpha),3),dtype = np.float)
    Read_cen_hal_MU = np.zeros((len(alpha),3),dtype = np.float)
    Read_cen_gla_GX = np.zeros((len(alpha),3),dtype = np.float)
    Read_cen_hal_GX = np.zeros((len(alpha),3),dtype = np.float)
    # save the data to figure the comparation fig
    m_s = np.zeros((T_scale,Nsize),dtype = np.float)
    m_g = np.zeros((T_scale,Nsize),dtype = np.float)
    m_s_gx = np.zeros((T_scale,Nsize),dtype = np.float)
    m_g_gx = np.zeros((T_scale,Nsize),dtype = np.float)
    rho_s = np.zeros((T_scale,Nsize-1),dtype = np.float)
    rho_g = np.zeros((T_scale,Nsize-1),dtype = np.float)
    rho_s_gx = np.zeros((T_scale,Nsize-1),dtype = np.float)
    rho_g_gx = np.zeros((T_scale,Nsize-1),dtype = np.float)
    m_rho_s = np.zeros((T_scale,Nsize),dtype = np.float)
    m_rho_g = np.zeros((T_scale,Nsize),dtype = np.float)
    m_rho_s_gx = np.zeros((T_scale,Nsize),dtype = np.float)
    m_rho_g_gx = np.zeros((T_scale,Nsize),dtype = np.float)
    z_mu = np.zeros(T_scale,dtype = np.float)
    z_gx = np.zeros(T_scale,dtype = np.float)
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
    id_z = com_z[iv[_id_]]
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
        snap_N = pygdr.readheader(
            '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
        print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f' % (
            snap_N[0], snap_N[1], snap_N[2], snap_N[3], snap_N[4], snap_N[5]))
        #snap_name = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap, 'header')
        # for type
        if snap_z <= id_z:
            z_mu[k] = snap_z
            snap_shot_gas = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_mass_gas = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
            snap_shot_DM = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
            try:
                snap_shot_star = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
                snap_mass_star = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
            except SystemExit:print('no star particles now')
            # for respective position
            snap_shot_bulge = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bulge')
            snap_shot_disk = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'disk')
            try:
                snap_shot_bndry = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bndry')
                snap_mass_BH = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','bndry')
            except SystemExit:
                print('no boundary particles now')
            # for density profile of gas(there only density about gas in the simulation data)
            # snap_dens = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'rho','gas')
            # the total mass distribution of snap_shot_128
            read_z = rnr.get_round_number_float(snap_z,3)
            # next,try to find the particles belong to halo 128000000000001(in MUSIC simulation)
            main_halo = asc.read(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%read_z,
                    converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z = changds.inv_chidas('%.3f'%read_z)
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
            # try to find the most densitive area in the central area(which will set as BCG,and the selection 
            # area as a sphere)
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
                compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+(cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
                ismin = find.find1d(compare_d,np.min(compare_d))
                cen_x_star = cen_po_star[ismin,0]
                cen_y_star = cen_po_star[ismin,1]
                cen_z_star = cen_po_star[ismin,2]  
                Read_cen_gla_MU[k,:] = np.array([cen_x_star,cen_y_star,cen_z_star])
                Read_cen_hal_MU[k,:] = np.array([x0,y0,z0])
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
                m_g[k,:] = inl_m_g
                if snap_N[-2] !=0:
                    inl_m_s = np.zeros(n_r,dtype = np.float)
                    r_star = np.sqrt((inlstar[:,0]-cen_x_star)**2+(inlstar[:,1]-cen_y_star)**2+
                                 (inlstar[:,2]-cen_z_star)**2)
                    for p in range(n_r):
                        ib = r_star <= r_bcg[p]
                        inl_m_s[p] = np.sum(inlmass_star[ib]*10**10)
                    m_s[k,:] = inl_m_s
                if snap_N[-1] !=0:
                    inl_BH = {}
                    inl_BH_m = np.zeros(n_r,dtype = np.float)
                    r_BH = np.sqrt((inlbndry[:,0]-cen_x_star)**2+(inlbndry[:,1]-cen_y_star)**2+
                                 (inlbndry[:,2]-cen_z_star)**2)
                    for p in range(n_r):
                        ic = r_BH <=r_bcg[p]
                        inl_BH[p] = inlbndry[ic,:]
                        inl_BH_m[p] = np.sum(inlmass_BH[ic]*10**10)
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
                rho_s[k,:] = inl_rho_s
                rho_g[k,:] = inl_rho_g
                inl_rho_b = inl_rho_g + inl_rho_s
                mean_rho_g = inl_m_g/(4*np.pi*r_bcg**3/3)
                m_rho_g[k,:] = mean_rho_g
                if snap_N[-2] !=0:
                    mean_rho_s = inl_m_s/(4*np.pi*r_bcg**3/3)
                    m_rho_s[k,:] = mean_rho_s
                mean_rho_b = inl_m_b/(4*np.pi*r_bcg**3/3)
                plt.figure()
                if snap_N[-2] != 0:
                    plt.hist2d(inlstar[:, 0], inlstar[:, 1], bins=[1000, 1000],
                       cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                    plt.scatter(cen_x_star, cen_y_star, s=20, c='m', marker='o',)
                    plt.scatter(x0, y0, s=20, c='k', marker='+')
                    plt.xlim(x0-R0, x0+R0)
                    plt.ylim(y0-R0, y0+R0)
                    plt.xticks([x0-R0,x0,x0+R0])
                    plt.yticks([y0, y0+R0], rotation=90)
                    plt.xlabel(r'$r-kpc/h$')
                    plt.ylabel(r'$r-kpc/h$')
                    plt.savefig(
                            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_star_mu_No_%.0f_z_%.0f.png' % (_id_, dd), dpi=600)
                    plt.show()
                    plt.close()
                plt.figure()
                # plt.hist2d(inlDM[:,0],inlDM[:,1],bins = [1000,1000],
                #           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
                # plt.hist2d(inlgas[:,0],inlgas[:,1],bins = [1000,1000],
                #           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    # plt.hist2d(inlstar[:,0],inlstar[:,1],bins = [1000,1000],
                    #       cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
                    plt.hist2d(inl_star[:, 0], inl_star[:, 1], bins=[num_bins, num_bins],
                       cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                # hsc.circles(cen_po_gas[0],cen_po_gas[1],s = resolution, c = 'k')
                # hsc.circles(cen_po_star[0],cen_po_star[1],s = resolution, c = 'm')
                # hsc.circles(x0,y0,s = 5*Nr*resolution ,c = 'b',alpha = 0.35)
                A = [[cen_x_star-10*Nr*resolution, cen_x_star +
                    10*Nr*resolution], [cen_x_star, cen_x_star]]
                A = np.array(A)
                B = [[cen_y_star, cen_y_star], [cen_y_star -
                    10*Nr*resolution, cen_y_star+10*Nr*resolution]]
                B = np.array(B)
                plt.plot(A[0], B[0], 'k-', lw=1, alpha=0.35,
                         label=r'$20kpc/h_{per line}-C_\ast$')
                plt.plot(A[1], B[1], 'k-', lw=1, alpha=0.35)
                C = [[x0-10*Nr*resolution, x0+10*Nr*resolution], [x0, x0]]
                D = [[y0, y0], [y0-10*Nr*resolution, y0+10*Nr*resolution]]
                plt.plot(C[0], D[0], 'b-', lw=2, alpha=0.35,
                         label=r'$20kpc/h_{per line}-C_h$')
                plt.plot(C[1], D[1], 'b-', lw=2, alpha=0.35)
                if snap_N[-1] !=0:
                    plt.scatter(inl_BH[n_r-1][:,0],inl_BH[n_r-1][:,1],s = 20,c = 'c',alpha = 0.35)
                plt.xlim(x0-R_range, x0+R_range)
                plt.ylim(y0-R_range, y0+R_range)
                plt.legend(loc=4, fontsize=5)
                plt.xlabel(r'$r-kpc/h$')
                plt.ylabel(r'$r-kpc/h$')
                plt.xticks([x0-R_range, x0, x0+R_range], size=10)
                plt.yticks([y0, y0+R_range], rotation=90, size=10)
                plt.title(r'$MUSIC_{%.3f}-resolution_{%.3f}$' % (snap_z, resolution))
                plt.savefig(
                        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_center_mu No_%.0f_z_%.0f.png' % (_id_, dd), dpi=600)
                plt.show()
                plt.close()
            except ValueError:
                print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
                break
    Read_cen_gla_MU = Read_cen_gla_MU[Read_cen_gla_MU != 0]
    Read_cen_gla_MU = Read_cen_gla_MU.reshape(np.int0(len(Read_cen_gla_MU)/3),3)
    Read_cen_hal_MU = Read_cen_hal_MU[Read_cen_hal_MU != 0]
    Read_cen_hal_MU = Read_cen_hal_MU.reshape(np.int0(len(Read_cen_hal_MU)/3),3) 
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_gla_%.0f_MU_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(Read_cen_gla_MU)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_gla_%.0f_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(len(Read_cen_gla_MU)):
            f['a'][t,:] = np.array(Read_cen_gla_MU[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_hal_MU.h5','w') as f:
        f['a'] = np.array(Read_cen_hal_MU)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_hal_MU.h5') as f:
        for t in range(len(Read_cen_hal_MU)):
            f['a'][t,:] = np.array(Read_cen_hal_MU[t,:])
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
        # section 3: read the data of GX simution
        snap_z_gx = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
        snap_z_gx = np.abs(snap_z_gx)
        print('Now redshift is %.3f' % snap_z_gx)
        snap_N_gx = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
        print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
    	      % (snap_N_gx[0], snap_N_gx[1], snap_N_gx[2], snap_N_gx[3], snap_N_gx[4], snap_N_gx[5]))
        #snap_name_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'header')
        # for type
        if snap_z_gx <= id_z_gx:
            z_gx[k] = snap_z_gx
            snap_shot_gas_gx = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_mass_gas_gx = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
            snap_shot_DM_gx = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
            try:
                snap_shot_star_gx = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
                snap_mass_star_gx = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
            except SystemExit:
                print('no star particles now')
            # for respective position
            snap_shot_bulge_gx = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bulge')
            snap_shot_disk_gx = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'disk')
            try:
                snap_shot_bndry_gx = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bndry')
                snap_mass_BH_gx = pygdr.readsnap(
                        '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','bndry')
            except SystemExit:
                print('no boundary particles now')
            # for density profile of gas(there only density about gas in the simulation data)
            # snap_dens_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'rho','gas')
            # the total mass distribution of snap_shot_128
            read_z_gx = rnr.get_round_number_float(snap_z_gx,3)
            # next,try to find the particles belong to halo 128000000000001(in GX simulation)
            main_halo_gx = asc.read(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%read_z_gx,
                    converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z_gx = changds.inv_chidas('%.3f'%read_z_gx)
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
                Read_cen_gla_GX[k,:] = np.array([cen_x_star_gx,cen_y_star_gx,cen_z_star_gx])
                Read_cen_hal_GX[k,:] = np.array([x0_gx,y0_gx,z0_gx])
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
                    inl_m_g_gx[q] = np.sum(inlmass_gas_gx[ia_gx]*10**10)
                m_g_gx[k,:] = inl_m_g_gx
                if snap_N_gx[-2] !=0:
                    inl_m_s_gx = np.zeros(n_r_gx,dtype = np.float)
                    r_star_gx = np.sqrt((inlstar_gx[:,0]-cen_x_star_gx)**2+(inlstar_gx[:,1]-cen_y_star_gx)**2+
                                 (inlstar_gx[:,2]-cen_z_star_gx)**2)
                    for q in range(n_r_gx):
                        ib_gx = r_star_gx <= r_bcg_gx[q]
                        inl_m_s_gx[q] = np.sum(inlmass_star_gx[ib_gx]*10**10)
                    m_s_gx[k,:] = inl_m_s_gx
                if snap_N_gx[-1] !=0:
                    inl_BH_gx = {}
                    inl_BH_m_gx = np.zeros(n_r_gx,dtype = np.float)
                    r_BH_gx = np.sqrt((inlbndry_gx[:,0]-cen_x_star_gx)**2+(inlbndry_gx[:,1]-cen_y_star_gx)**2+
                                 (inlbndry_gx[:,2]-cen_z_star_gx)**2)
                    for q in range(n_r_gx):
                        ic_gx = r_BH_gx <=r_bcg_gx[q]
                        inl_BH_gx[q] = inlbndry_gx[ic_gx,:]
                        inl_BH_m_gx[q] = np.sum(inlmass_BH_gx[ic_gx]*10**10)
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
                rho_s_gx[k,:] = inl_rho_s_gx
                rho_g_gx[k,:] = inl_rho_g_gx
                inl_rho_b_gx = inl_rho_g_gx + inl_rho_s_gx
                mean_rho_g_gx = inl_m_g_gx/(4*np.pi*r_bcg_gx**3/3)
                m_rho_g_gx[k,:] = mean_rho_g_gx
                if snap_N_gx[-2] !=0:
                    mean_rho_s_gx = inl_m_s_gx/(4*np.pi*r_bcg_gx**3/3)
                    m_rho_s_gx[k,:] = mean_rho_g_gx
                mean_rho_b_gx = inl_m_b_gx/(4*np.pi*r_bcg_gx**3/3)
                plt.figure()
                if snap_N_gx[-2] != 0:
                    plt.hist2d(inlstar_gx[:, 0], inlstar_gx[:, 1], bins=[1000, 1000],
                          cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                    plt.scatter(cen_x_star_gx,cen_y_star_gx, s=20, c='m', marker='o',)
                    plt.scatter(x0_gx, y0_gx, s=20, c='k', marker='+')
                    plt.xlim(x0_gx-R0_gx, x0_gx+R0_gx)
                    plt.ylim(y0_gx-R0_gx, y0_gx+R0_gx)
                    plt.xticks([x0_gx-R0_gx,x0_gx,x0_gx+R0_gx])
                    plt.yticks([y0_gx, y0_gx+R0_gx], rotation=90)
                    plt.xlabel(r'$r-kpc/h$')
                    plt.ylabel(r'$r-kpc/h$')
                    plt.savefig(
                            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_star_gx_No_%.0f_z_%.0f.png' % (_id_, dd), dpi=600)
                    plt.show()
                    plt.close()
                plt.figure()
                # plt.hist2d(inlDM_gx[:,0],inlDM_gx[:,1],bins = [1000,1000],
                #           cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
                # plt.hist2d(inlgas_gx[:,0],inlgas_gx[:,1],bins = [1000,1000],
                #           cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
                if snap_N_gx[-2] != 0:
                   # plt.hist2d(inlstar_gx[:,0],inlstar_gx[:,1],bins = [1000,1000],
                   #       cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
                   plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 1], bins=[num_bins_gx, num_bins_gx],
                           cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
                   plt.colorbar(label=r'$log[Number density]$')
                # hsc.circles(cen_po_gas_gx[0],cen_po_gas_gx[1],s = resolution, c = 'k')
                # hsc.circles(cen_po_star_gx[0],cen_po_star_gx[1],s = resolution, c = 'm')
                # hsc.circles(x0_gx,y0_gx,s = 5*Nr*resolution,c = 'b',alpha = 0.35)
                A_gx = [[cen_x_star_gx-10*Nr*resolution, cen_x_star_gx +
                    10*Nr*resolution], [cen_x_star_gx, cen_x_star_gx]]
                A_gx = np.array(A_gx)
                B_gx = [[cen_y_star_gx, cen_y_star_gx], 
                        [cen_y_star_gx-10*Nr*resolution, cen_y_star_gx+10*Nr*resolution]]
                B_gx = np.array(B_gx)
                plt.plot(A_gx[0], B_gx[0], 'k-', lw=1, alpha=0.35,
                         label=r'$20kpc/h_{per line}-C_\ast$')
                plt.plot(A_gx[1], B_gx[1], 'k-', lw=1, alpha=0.35)
                C_gx = [[x0_gx-10*Nr*resolution, x0_gx+10*Nr*resolution], [x0_gx, x0_gx]]
                D_gx = [[y0_gx, y0_gx], [y0_gx-10*Nr*resolution, y0_gx+10*Nr*resolution]]
                plt.plot(C_gx[0], D_gx[0], 'b-', lw=2, alpha=0.35,
                         label=r'$20kpc/h_{per line}-C_h$')
                plt.plot(C_gx[1], D_gx[1], 'b-', lw=2, alpha=0.35)
                if snap_N_gx[-1] !=0:
                    plt.scatter(inl_BH_gx[n_r_gx-1][:,0],inl_BH_gx[n_r_gx-1][:,1],s = 20,c = 'c',alpha = 0.35)
                plt.xlim(x0_gx-R_range_gx, x0_gx+R_range_gx)
                plt.ylim(y0_gx-R_range_gx, y0_gx+R_range_gx)
                plt.legend(loc=4, fontsize=5)
                plt.xlabel(r'$r-kpc/h$')
                plt.ylabel(r'$r-kpc/h$')
                plt.xticks([x0_gx-R_range_gx, x0_gx, x0_gx+R_range_gx], size=10)
                plt.yticks([y0_gx, y0_gx+R_range_gx], rotation=90, size=10)
                plt.title(r'$GadgetX_{%.3f}-resolution_{%.3f}$' %(snap_z_gx, resolution))
                plt.savefig(
                        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_Multi_center_gx_No_%.0f_z_%.0f.png' %(_id_, dd), dpi=600)
                plt.show()
                plt.close()
            except ValueError:
                print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z_gx)
                break
    Read_cen_gla_GX = Read_cen_gla_GX[Read_cen_gla_GX != 0]
    Read_cen_gla_GX = Read_cen_gla_GX.reshape(np.int0(len(Read_cen_gla_GX)/3),3)
    Read_cen_hal_GX = Read_cen_hal_GX[Read_cen_hal_GX != 0]
    Read_cen_hal_GX = Read_cen_hal_GX.reshape(np.int0(len(Read_cen_hal_GX)/3),3)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_gla_%.0f_GX_%.3f.h5'%(size_BCG,resolution),'w') as f:
        f['a'] = np.array(Read_cen_gla_GX)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_gla_%.0f_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        for t in range(len(Read_cen_gla_GX)):
            f['a'][t,:] = np.array(Read_cen_gla_GX[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_hal_GX.h5','w') as f:
        f['a'] = np.array(Read_cen_hal_GX)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_cen_hal_GX.h5') as f:
        for t in range(len(Read_cen_hal_GX)):
            f['a'][t,:] = np.array(Read_cen_hal_GX[t,:])
    il = m_s[:,0] != 0
    zmu = np.array(z_mu[il])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_z_mu.h5','w') as f:
        f['a'] = np.array(zmu)
    m_s = np.array(m_s[il,:])
    m_g = np.array(m_g[il,:])
    m_b = m_s + m_g
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_s_mu.h5','w') as f:
        f['a'] = np.array(m_s)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_s_mu.h5') as f:
        for t in range(np.int0(m_s.shape[0])):
            f['a'][t,:] = np.array(m_s[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_g_mu.h5','w') as f:
        f['a'] = np.array(m_g)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_g_mu.h5') as f:
        for t in range(np.int0(m_g.shape[0])):
            f['a'][t,:] = np.array(m_g[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_b_mu.h5','w') as f:
        f['a'] = np.array(m_b)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_b_mu.h5') as f:
        for t in range(np.int0(m_b.shape[0])):
            f['a'][t,:] = np.array(m_b[t,:])
    rho_s = np.array(rho_s[il,:])
    rho_g = np.array(rho_g[il,:])
    rho_b = rho_s + rho_g
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_s_mu.h5','w') as f:
        f['a'] = np.array(rho_s)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_s_mu.h5') as f:
        for t in range(np.int0(rho_s.shape[0])):
            f['a'][t,:] = np.array(rho_s[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_g_mu.h5','w') as f:
        f['a'] = np.array(rho_g)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_g_mu.h5') as f:
        for t in range(np.int0(rho_g.shape[0])):
            f['a'][t,:] = np.array(rho_g[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_b_mu.h5','w') as f:
        f['a'] = np.array(rho_b)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_b_mu.h5') as f:
        for t in range(np.int0(rho_b.shape[0])):
            f['a'][t,:] = np.array(rho_b[t,:])
    m_rho_s = np.array(m_rho_s[il,:])
    m_rho_g = np.array(m_rho_g[il,:])
    m_rho_b = m_rho_s + m_rho_g
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_s_mu.h5','w') as f:
        f['a'] = np.array(m_rho_s)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_s_mu.h5') as f:
        for t in range(np.int0(m_rho_s.shape[0])):
            f['a'][t,:] = np.array(m_rho_s[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_g_mu.h5','w') as f:
        f['a'] = np.array(m_rho_g)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_g_mu.h5') as f:
        for t in range(np.int0(m_rho_g.shape[0])):
            f['a'][t,:] = np.array(m_rho_g[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_b_mu.h5','w') as f:
        f['a'] = np.array(m_rho_b)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_b_mu.h5') as f:
        for t in range(np.int0(m_rho_b.shape[0])):
            f['a'][t,:] = np.array(m_rho_b[t,:])
    il_gx = m_s_gx[:,0] != 0
    zgx = np.array(z_gx[il_gx])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_z_gx.h5','w') as f:
        f['a'] = np.array(zgx)
    m_s_gx = np.array(m_s_gx[il_gx,:])
    m_g_gx = np.array(m_g_gx[il_gx,:])
    m_b_gx = m_s_gx + m_g_gx
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_s_gx.h5','w') as f:
        f['a'] = np.array(m_s_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_s_gx.h5') as f:
        for t in range(np.int0(m_s_gx.shape[0])):
            f['a'][t,:] = np.array(m_s_gx[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_g_gx.h5','w') as f:
        f['a'] = np.array(m_g_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_g_gx.h5') as f:
        for t in range(np.int0(m_g_gx.shape[0])):
            f['a'][t,:] = np.array(m_g_gx[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_b_gx.h5','w') as f:
        f['a'] = np.array(m_b_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_b_gx.h5') as f:
        for t in range(np.int0(m_b_gx.shape[0])):
            f['a'][t,:] = np.array(m_b[t,:])
    rho_s_gx = np.array(rho_s_gx[il_gx,:])
    rho_g_gx = np.array(rho_g_gx[il_gx,:])
    rho_b_gx = rho_s_gx + rho_g_gx
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_s_gx.h5','w') as f:
        f['a'] = np.array(rho_s_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_s_gx.h5') as f:
        for t in range(np.int0(rho_s_gx.shape[0])):
            f['a'][t,:] = np.array(rho_s_gx[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_g_gx.h5','w') as f:
        f['a'] = np.array(rho_g_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_g_gx.h5') as f:
        for t in range(np.int0(rho_g_gx.shape[0])):
            f['a'][t,:] = np.array(rho_g_gx[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_b_gx.h5','w') as f:
        f['a'] = np.array(rho_b_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_b_gx.h5') as f:
        for t in range(np.int0(rho_b_gx.shape[0])):
            f['a'][t,:] = np.array(rho_b_gx[t,:])
    m_rho_s_gx = np.array(m_rho_s_gx[il_gx,:])
    m_rho_g_gx = np.array(m_rho_g_gx[il_gx,:])
    m_rho_b_gx = m_rho_s_gx + m_rho_g_gx
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_s_gx.h5','w') as f:
        f['a'] = np.array(m_rho_s_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_s_gx.h5') as f:
        for t in range(np.int0(m_rho_s_gx.shape[0])):
            f['a'][t,:] = np.array(m_rho_s_gx[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_g_gx.h5','w') as f:
        f['a'] = np.array(m_rho_g_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_g_gx.h5') as f:
        for t in range(np.int0(m_rho_g_gx.shape[0])):
            f['a'][t,:] = np.array(m_rho_g_gx[t,:])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_b_gx.h5','w') as f:
        f['a'] = np.array(m_rho_b_gx)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_b_gx.h5') as f:
        for t in range(np.int0(m_rho_b_gx.shape[0])):
            f['a'][t,:] = np.array(m_rho_b_gx[t,:])
    return zgx
#find_BCG(resolution = True,g_size = True)
