# this file try to read the temperature distribution of gas for all of the snap_shot
import matplotlib as mpl
mpl.use('Agg')
from pygadgetreader import *
import h5py
import scipy.stats as st
import find
import pygadgetreader as pygdr
import numpy as np
import astropy.io.ascii as asc
import matplotlib.pyplot as plt
import pandas as pd
import astropy.constants as C
import astropy.units as U
import find
#chang the units
F = C.k_B
M_H = C.u.value
M = U.M_sun
L = U.km/U.s
Ms = C.M_sun.value
G = C.G
#change m,km to kpc,Mpc
C1_m = U.m.to(U.kpc)
C1_km = U.km.to(U.kpc)
C2_m = U.m.to(U.Mpc)
C2_km = U.km.to(U.Mpc)
W = G.to((U.km)**3*U.s**(-2)/(M))
xH = 0.76 #constant used to calculate the temperature
#handle the data
resolution = 1
size_BCG = 100
_id_ = 0
Nr = 1./resolution
# to make sure the bins density is the same
N_size = 50.*(size_BCG/100.0)
Nsize = np.int0(np.ceil(N_size))
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)
#read the data of megertree 
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
#next read and save the data of Temperature and SFR
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
            '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'redshift')
    snap_z = np.abs(snap_z)
    if snap_z <= id_z:
        print('Now redshift is %.3f'%snap_z)
        snap_N = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'npartTotal')
        print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
              %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
        snap_name = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'header')
        print(snap_name)
        
        snap_mass_gas = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
        snap_inE = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'U','gas') 
        snap_rho = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'RHO','gas') 
        snap_OmegaB = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'Omega0')
        print('omegaB=',snap_OmegaB)
        OmegaB = snap_OmegaB
        snap_H = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'h')
        snap_time = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'time')
        time = np.abs(snap_time)
        #PS:NE -- electron number friction is exactlly for each particle,and the same as NH.
        snap_NE = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'NE','gas')
        Ne = snap_NE
        snap_NH = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'NH','gas')
        Nh = snap_NH
        
        yhe = (1-xH)/(4*xH)
        mean_m_weight = (1+4*yhe)/(1+yhe+Ne)
        V_unit = 10**5*np.sqrt(time)
        K_B = F.value*10**7
        mp = M_H*10**3
        snap_T = snap_inE*(5/3-1)*V_unit**2*mp*mean_m_weight/K_B
        H = snap_H*100
        rho_bar = OmegaB*3*H**2/(8*np.pi*W.value)
        ratio_rho = snap_rho/rho_bar
        
        array1 = np.array([snap_T,ratio_rho])
        with h5py.File('ratio_T_rho_h%.0f_z%.3f_MU.h5'%(_id_,snap_z),'w') as f:
            f['a'] = array1
        with h5py.File('ratio_T_rho_h%.0f_z%.3f_MU.h5'%(_id_,snap_z)) as f:
            for la in range(array1.shape[0]):
                f['a'][la,:]= array1[la,:]
        
        #read the position of gas
        snap_shot_gas = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'pos','gas')
        try:
            snap_shot_star = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'pos','star')
            snap_mass_star = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
        except SystemExit:
            print('no star particles now')    
        snap_SFR = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'SFR ',0)
        array2 = np.array([snap_rho,snap_T,snap_SFR])
        with h5py.File('rho_T_h%.0f_z%.3f_MU.h5'%(_id_,snap_z),'w') as f:
            f['a'] = array2
        with h5py.File('rho_T_h%.0f_z%.3f_MU.h5'%(_id_,snap_z)) as f:
            for la in range(array2.shape[0]):
                f['a'][la,:] = array2[la,:]
                
        #the total mass distribution of snap_shot_128
        main_halo = asc.read(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z,
                             converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
        goal_z = np.float('%.3f'%snap_z)
        goalz = find.find1d(com_z,goal_z)
        goal_halo = tree_line[_id_,goalz]
        #check the halo in this redshift
        try:
            HALO = np.array(main_halo['col1'])
            check_id = find.find1d(HALO,goal_halo) 
            Rvir = np.array(main_halo['col12'])
            xhalo = np.array(main_halo['col6'])
            yhalo = np.array(main_halo['col7'])
            zhalo = np.array(main_halo['col8'])
            Nbins = np.array(main_halo['col37'])
            x0 = xhalo[check_id]
            y0 = yhalo[check_id]
            z0 = zhalo[check_id]
            R0 = Rvir[check_id]
            nbins0 = Nbins[check_id]
            #figure the temperature distribution of main halo 
            main_profile = pd.read_table(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%snap_z,
                                         dtype = np.float)
            r = np.array(main_profile['#r(1)'])
            r = np.abs(r)
            if check_id == 0:
                R = r[0:nbins0]
            else:
                sum_bin = np.sum(Nbins[0:check_id])
                R = r[sum_bin : sum_bin + nbins0]
            L_R = len(R)
            dgas = np.sqrt((snap_shot_gas[:, 0]-x0)**2 +
                           (snap_shot_gas[:, 1]-y0)**2 + (snap_shot_gas[:, 2]-z0)**2)
            ig = dgas <= R0
            inlgas = snap_shot_gas[ig, :]
            inlmass_gas = snap_mass_gas[ig]
            
            inl_rho = snap_rho[ig]
            inl_inE = snap_inE[ig]
            inl_T = snap_T[ig]
            halo_SFR = snap_SFR[ig]
            array3 = np.array([inl_rho,halo_SFR,inl_T])
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/main%.0f_halo_rho_T_z%.3f_MU.h5'%(_id_,snap_z),'w') as f:
                f['a'] = array3
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/main%.0f_halo_rho_T_z%.3f_MU.h5'%(_id_,snap_z)) as f:
                for la in range(array3.shape[0]):
                    f['a'][la,:] = array3[la,:]
                    
            #to get the mean temperature
            mean_T_m = np.zeros(L_R,dtype = np.float)
            mean_T_n = np.zeros(L_R,dtype = np.float)
            for p in range(L_R):
                if p == 0:
                    igk = dgas <= R[p]
                    gasmass = snap_mass_gas[igk]
                    if len(gasmass) == 0:
                        mean_T_m[p] = 0.
                        mean_T_n[p] = 0.
                    else:
                        gasT = snap_T[igk]
                        mean_T_m[p] = np.sum(gasT*gasmass*10**10*Ms)/(np.sum(gasmass*10**10*Ms))
                        mean_T_n[p] = np.sum(gasT)/np.float(len(gasmass))
                else:
                    igk = (dgas <= R[p]) & (dgas > R[p-1])
                    gasmass = snap_mass_gas[igk]
                    if len(gasmass) == 0:
                        mean_T_m[p] = 0.
                        mean_T_n[p] = 0.
                    else:
                        gasT = snap_T[igk]
                        mean_T_m[p] = np.sum(gasT*gasmass*10**10*Ms)/(np.sum(gasmass*10**10*Ms))
                        mean_T_n[p] = np.sum(gasT)/np.float(len(gasmass))
            array4 = np.array([R,mean_T_m,mean_T_n])
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/meanT_R_h%.0f_z%.3f_MU.h5'%(_id_,snap_z),'w') as f:
                f['a'] = array4
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/meanT_R_h%.0f_z%.3f_MU.h5'%(_id_,snap_z)) as f:
                for la in range(array4.shape[0]):
                    f['a'][la,:] = array4[la,:]
            
            if snap_N[-2] != 0:
                dstar = np.sqrt((snap_shot_star[:, 0]-x0)**2+(snap_shot_star[:, 1]-y0)**2 + (snap_shot_star[:, 2]-z0)**2)
                ids = dstar <= R0
                inlstar = snap_shot_star[ids, :]
                inlmass_star = snap_mass_star[ids]
            R_range = 150.0
            edge_x = np.array([x0-R_range,x0+R_range])
            edge_y = np.array([y0-R_range,y0+R_range])
            edge_z = np.array([z0-R_range,z0+R_range])
            iv_s = ((inlstar[:,0]<=edge_x[1])&(inlstar[:,0]>=edge_x[0]))&(
                    (inlstar[:,1]<=edge_y[1])&(inlstar[:,1]>=edge_y[0]))&((inlstar[:,2]<=edge_z[1])&(inlstar[:,2]>=edge_z[0]))
            test_star = inlstar[iv_s,:]
            iv_g = ((inlgas[:,0]<=edge_x[1])&(inlgas[:,0]>=edge_x[0]))&(
                    (inlgas[:,1]<=edge_y[1])&(inlgas[:,1]>=edge_y[0]))&((inlgas[:,2]<=edge_z[1])&(inlgas[:,2]>=edge_z[0]))
            test_gas = inlgas[iv_g,:]
            # central postion and density of star
            num_bins = np.ceil(R_range*2/resolution)
            try:
                hist_star, edge_star = np.histogramdd(test_star, bins=(num_bins, num_bins, num_bins))
                bin_x_star = np.array(edge_star[0])
                bin_y_star = np.array(edge_star[1])
                bin_z_star = np.array(edge_star[2])
                inumber1 = hist_star >= 10
                #to find the first five density center and then choose the closed one to the halo center
                maxN = len(inumber1)
                cen_po_star = np.zeros((maxN,3),dtype = np.float)
                hist_use = hist_star
                for q in range(maxN):
                    is_max = np.unravel_index(np.argmax(hist_use, axis=None), hist_use.shape)
                    cenxstar = (bin_x_star[is_max[0]+1]+bin_x_star[is_max[0]])/2.0
                    cenystar = (bin_y_star[is_max[1]+1]+bin_y_star[is_max[1]])/2.0
                    cenzstar = (bin_z_star[is_max[2]+1]+bin_z_star[is_max[2]])/2.0
                    cen_po_star[q,:] = np.array([cenxstar,cenystar,cenzstar])
                    hist_use[is_max] = 0.0
                compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+(cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
                ismin = find.find1d(compare_d,np.min(compare_d))
                cen_x_star = cen_po_star[ismin,0]
                cen_y_star = cen_po_star[ismin,1]
                cen_z_star = cen_po_star[ismin,2]
                dgas_galaxy = np.sqrt((snap_shot_gas[:, 0]-cen_x_star)**2 +
                              (snap_shot_gas[:, 1]-cen_y_star)**2 + (snap_shot_gas[:, 2]-cen_z_star)**2) 
                i_galaxy = dgas_galaxy <= size_BCG
                T_gas = snap_T[i_galaxy]
                rho_gas = snap_rho[i_galaxy]
                ## figure the temperature profile of BCG
                R_BCG = size_BCG
                r_bcg = np.logspace(1e-3,np.log10(R_BCG),Nsize)
                n_r = len(r_bcg)
                T_galaxy_m = np.zeros(n_r,dtype = np.float)
                T_galaxy_n = np.zeros(n_r,dtype = np.float)
                for q in range(n_r):
                    if q == 0:
                        igp = dgas_galaxy <= r_bcg[q]
                        gasm = snap_mass_gas[igp]
                        if len(gasm) == 0:
                            T_galaxy_m[q] = 0.
                            T_galaxy_n[q] = 0.
                        else:
                            gasTg = snap_T[igp]
                            T_galaxy_m[q] = np.sum(gasTg*gasm*10**10*Ms)/(np.sum(gasm*10**10*Ms))
                            T_galaxy_n[q] = np.sum(gasTg)/np.float(len(gasm))
                    else:
                        igp = (dgas_galaxy <= r_bcg[q]) & (dgas_galaxy > r_bcg[q-1])
                        gasm = snap_mass_gas[igp]
                        if len(gasm) == 0:
                            T_galaxy_m[q] = 0.
                            T_galaxy_n[q] = 0.
                        else:
                            gasTg = snap_T[igp]
                            T_galaxy_m[q] = np.sum(gasTg*gasm*10**10*Ms)/(np.sum(gasm*10**10*Ms))
                            T_galaxy_n[q] = np.sum(gasTg)/np.float(len(gasm))
                T_galaxy_n[np.isnan(T_galaxy_n)] = 0.
                T_galaxy_m[np.isnan(T_galaxy_m)] = 0.
                array5 = np.array([r_bcg,T_galaxy_m,T_galaxy_n])
                with h5py.File(
                        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_T_R_h%.0f_z%.3f_MU.h5'%(_id_,snap_z),'w') as f:
                    f['a'] = array5
                with h5py.File(
                        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_T_R_h%.0f_z%.3f_MU.h5'%(_id_,snap_z)) as f:
                    for la in range(array5.shape[0]):
                        f['a'][la,:] = array5[la,:]
            except ValueError:
                print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
        except IndexError:
            print('Now is no right merger tree to trace back!') 
#for the gadgetX simulation
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
            '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'redshift')
    snap_z = np.abs(snap_z)
    if snap_z <= id_z_gx:
        print('Now redshift is %.3f'%snap_z)
        snap_N = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'npartTotal')
        print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
              %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
        snap_name = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'header')
        print(snap_name)
        
        snap_mass_gas = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','gas')
        snap_inE = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'U','gas') 
        snap_rho = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'RHO','gas') 
        snap_OmegaB = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'Omega0')
        print('omegaB=',snap_OmegaB)
        OmegaB = snap_OmegaB
        snap_H = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'h')
        snap_time = pygdr.readheader(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'time')
        time = np.abs(snap_time)
        #PS:NE -- electron number friction is exactlly for each particle,and the same as NH.
        snap_NE = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'NE','gas')
        Ne = snap_NE
        snap_NH = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'NH','gas')
        Nh = snap_NH
        ### try to get the temperature from U
        yhe = (1-xH)/(4*xH)
        mean_m_weight = (1+4*yhe)/(1+yhe+Ne)
        V_unit = 10**5*np.sqrt(time)
        K_B = F.value*10**7
        mp = M_H*10**3
        snap_T_confer = snap_inE*(5/3-1)*V_unit**2*mp*mean_m_weight/K_B 
        H = snap_H*100
        rho_bar = OmegaB*3*H**2/(8*np.pi*W.value)
        ratio_rho = snap_rho/rho_bar
        
        #read the gas position and temperature
        try:
            snap_T = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'TEMP','gas')
        except SystemExit:
            print('no gas temperature now')
        #read the position of gas
        snap_shot_gas = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'pos','gas')
        try:
            snap_shot_star = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'pos','star')
            snap_mass_star = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
        except SystemExit:
            print('no star particles now')    
        snap_SFR = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'SFR ',0)
        array1_gx = np.array([snap_T,ratio_rho])
        with h5py.File('ratio_T_rho_h%.0f_z%.3f_GX.h5'%(_id_,snap_z),'w') as f:
            f['a'] = array1_gx
        with h5py.File('ratio_T_rho_h%.0f_z%.3f_GX.h5'%(_id_,snap_z)) as f:
            for la in range(array1_gx.shape[0]):
                f['a'][la,:]= array1_gx[la,:]
        
        array2_gx = np.array([snap_rho,snap_T,snap_SFR])
        with h5py.File('rho_T_h%.0f_z%.3f_GX.h5'%(_id_,snap_z),'w') as f:
            f['a'] = array2_gx
        with h5py.File('rho_T_h%.0f_z%.3f_GX.h5'%(_id_,snap_z)) as f:
            for la in range(array2_gx.shape[0]):
                f['a'][la,:] = array2_gx[la,:]
        
        #the total mass distribution of snap_shot_128
        main_halo_gx = asc.read(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z,
                             converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
        goal_z_gx = np.float('%.3f'%snap_z)
        goalz_gx = find.find1d(com_z_gx,goal_z_gx)
        goal_halo_gx = tree_line_gx[_id_,goalz_gx]
        #check the halo in this redshift
        try:
            HALO_gx = np.array(main_halo_gx['col1'])
            check_id = find.find1d(HALO_gx,goal_halo_gx) 
            Rvir = np.array(main_halo_gx['col12'])
            xhalo = np.array(main_halo_gx['col6'])
            yhalo = np.array(main_halo_gx['col7'])
            zhalo = np.array(main_halo_gx['col8'])
            Nbins = np.array(main_halo_gx['col37'])
            x0 = xhalo[check_id]
            y0 = yhalo[check_id]
            z0 = zhalo[check_id]
            R0 = Rvir[check_id]
            nbins0 = Nbins[check_id]
            #figure the temperature distribution of main halo 
            main_profile = pd.read_table(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%snap_z,
                                         dtype = np.float)
            r = np.array(main_profile['#r(1)'])
            r = np.abs(r)
            if check_id == 0:
                R = r[0:nbins0]
            else:
                sum_bin = np.sum(Nbins[0:check_id])
                R = r[sum_bin : sum_bin + nbins0]
            L_R = len(R)
            dgas = np.sqrt((snap_shot_gas[:, 0]-x0)**2 +
                           (snap_shot_gas[:, 1]-y0)**2 + (snap_shot_gas[:, 2]-z0)**2)
            ig = dgas <= R0
            inlgas = snap_shot_gas[ig, :]
            inlmass_gas = snap_mass_gas[ig]
            
            inl_rho = snap_rho[ig]
            inl_inE = snap_inE[ig]
            inl_T = snap_T[ig]
            halo_SFR = snap_SFR[ig]
    
            array3_gx = np.array([inl_rho,halo_SFR,inl_T])
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/main%.0f_halo_rho_T_z%.3f_GX.h5'%(_id_,snap_z),'w') as f:
                f['a'] = array3_gx
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/main%.0f_halo_rho_T_z%.3f_GX.h5'%(_id_,snap_z)) as f:
                for la in range(array3_gx.shape[0]):
                    f['a'][la,:] = array3_gx[la,:]
            #to get the mean temperature
            mean_T_m = np.zeros(L_R,dtype = np.float)
            mean_T_n = np.zeros(L_R,dtype = np.float)
            for p in range(L_R):
                if p == 0:
                    igk = dgas <= R[p]
                    gasmass = snap_mass_gas[igk]
                    if len(gasmass) == 0:
                        mean_T_m[p] = 0.
                        mean_T_n[p] = 0.
                    else:
                        gasT = snap_T[igk]
                        mean_T_m[p] = np.sum(gasT*gasmass*10**10*Ms)/(np.sum(gasmass*10**10*Ms))
                        mean_T_n[p] = np.sum(gasT)/np.float(len(gasmass))
                else:
                    igk = (dgas <= R[p]) & (dgas > R[p-1])
                    gasmass = snap_mass_gas[igk]
                    if len(gasmass) == 0:
                        mean_T_m[p] = 0.
                        mean_T_n[p] = 0.
                    else:
                        gasT = snap_T[igk]
                        mean_T_m[p] = np.sum(gasT*gasmass*10**10*Ms)/(np.sum(gasmass*10**10*Ms))
                        mean_T_n[p] = np.sum(gasT)/np.float(len(gasmass))
            array4_gx = np.array([R,mean_T_m,mean_T_n])
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/meanT_R_h%.0f_z%.3f_GX.h5'%(_id_,snap_z),'w') as f:
                f['a'] = array4_gx
            with h5py.File(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/meanT_R_h%.0f_z%.3f_GX.h5'%(_id_,snap_z)) as f:
                for la in range(array4_gx.shape[0]):
                    f['a'][la,:] = array4_gx[la,:]
            
            if snap_N[-2] != 0:
                dstar = np.sqrt((snap_shot_star[:, 0]-x0)**2+(snap_shot_star[:, 1]-y0)**2 + (snap_shot_star[:, 2]-z0)**2)
                ids = dstar <= R0
                inlstar = snap_shot_star[ids, :]
                inlmass_star = snap_mass_star[ids]
            R_range = 150.0
            edge_x = np.array([x0-R_range,x0+R_range])
            edge_y = np.array([y0-R_range,y0+R_range])
            edge_z = np.array([z0-R_range,z0+R_range])
            iv_s = ((inlstar[:,0]<=edge_x[1])&(inlstar[:,0]>=edge_x[0]))&(
                    (inlstar[:,1]<=edge_y[1])&(inlstar[:,1]>=edge_y[0]))&((inlstar[:,2]<=edge_z[1])&(inlstar[:,2]>=edge_z[0]))
            test_star = inlstar[iv_s,:]
            iv_g = ((inlgas[:,0]<=edge_x[1])&(inlgas[:,0]>=edge_x[0]))&(
                    (inlgas[:,1]<=edge_y[1])&(inlgas[:,1]>=edge_y[0]))&((inlgas[:,2]<=edge_z[1])&(inlgas[:,2]>=edge_z[0]))
            test_gas = inlgas[iv_g,:]
            # central postion and density of star
            num_bins = np.ceil(R_range*2/resolution)
            try:
                hist_star, edge_star = np.histogramdd(test_star, bins=(num_bins, num_bins, num_bins))
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
                dgas_galaxy = np.sqrt((snap_shot_gas[:, 0]-cen_x_star)**2 +
                              (snap_shot_gas[:, 1]-cen_y_star)**2 + (snap_shot_gas[:, 2]-cen_z_star)**2) 
                i_galaxy = dgas_galaxy <= size_BCG
                T_gas = snap_T[i_galaxy]
                rho_gas = snap_rho[i_galaxy]
                ## figure the temperature profile of BCG
                R_BCG = size_BCG
                r_bcg = np.logspace(1e-3,np.log10(R_BCG),Nsize)
                n_r = len(r_bcg)
                T_galaxy_m = np.zeros(n_r,dtype = np.float)
                T_galaxy_n = np.zeros(n_r,dtype = np.float)
                for q in range(n_r):
                    if q == 0:
                        igp = dgas_galaxy <= r_bcg[q]
                        gasm = snap_mass_gas[igp]
                        if len(gasm) == 0:
                            T_galaxy_m[q] = 0.
                            T_galaxy_n[q] = 0.
                        else:
                            gasTg = snap_T[igp]
                            T_galaxy_m[q] = np.sum(gasTg*gasm*10**10*Ms)/(np.sum(gasm*10**10*Ms))
                            T_galaxy_n[q] = np.sum(gasTg)/np.float(len(gasm))
                    else:
                        igp = (dgas_galaxy <= r_bcg[q]) & (dgas_galaxy > r_bcg[q-1])
                        gasm = snap_mass_gas[igp]
                        if len(gasm) == 0:
                            T_galaxy_m[q] = 0.
                            T_galaxy_n[q] = 0.
                        else:
                            gasTg = snap_T[igp]
                            T_galaxy_m[q] = np.sum(gasTg*gasm*10**10*Ms)/(np.sum(gasm*10**10*Ms))
                            T_galaxy_n[q] = np.sum(gasTg)/np.float(len(gasm))
                T_galaxy_n[np.isnan(T_galaxy_n)] = 0.
                T_galaxy_m[np.isnan(T_galaxy_m)] = 0.
                array5_gx = np.array([r_bcg,T_galaxy_m,T_galaxy_n])
                with h5py.File(
                        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_T_R_h%.0f_z%.3f_GX.h5'%(_id_,snap_z),'w') as f:
                    f['a'] = array5_gx
                with h5py.File(
                        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_T_R_h%.0f_z%.3f_GX.h5'%(_id_,snap_z)) as f:
                    for la in range(array5_gx.shape[0]):
                        f['a'][la,:] = array5_gx[la,:]            
            except ValueError:
                print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
        except IndexError:
            print('Now is no right merger tree to trace back!')
