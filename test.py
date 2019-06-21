# this file try to find the BCG in a given halo
import matplotlib as mpl
#mpl.use('Agg')
#from pygadgetreader import *
#import h5py
#import find
import pygadgetreader as pygdr
import numpy as np
import astropy.io.ascii as asc
import matplotlib.pyplot as plt
import pandas as pd
import astropy.constants as C
import astropy.units as U
import handy
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
No_snap = '128'
_id_ = 0
resolution = 1.0
Nr = 1./resolution
size_BCG = 100
#to make sure the bins density is the same
N_size = 50.*(size_BCG/100.0)
Nsize = np.ceil(N_size)
choose = 1
#######################section1:read out the position data of MUSIC
if choose == 0:
    snap_z = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'redshift')
    snap_z = np.abs(snap_z)
    print('Now redshift is %.3f'%snap_z)
    snap_N = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'npartTotal')
    print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
          %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
    snap_name = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'header')
    print(snap_name)
    
    snap_mass_gas = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'mass','gas')
    snap_inE = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'U','gas') 
    snap_rho = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'RHO','gas') 
    snap_OmegaB = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'Omega0')
    print('omegaB=',snap_OmegaB)
    OmegaB = snap_OmegaB
    snap_H = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'h')
    snap_time = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'time')
    time = np.abs(snap_time)
    #PS:NE -- electron number friction is exactlly for each particle,and the same as NH.
    snap_NE = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'NE','gas')
    Ne = snap_NE
    snap_NH = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'NH','gas')
    Nh = snap_NH
    '''
    yhe = (1-xH)/(4*xH)
    mean_m_weight = (1+4*yhe)/(1+yhe+Ne)
    V_unit = 10**5*np.sqrt(time)
    K_B = F.value*10**7
    mp = M_H*10**3
    snap_T = snap_inE*(5/3-1)*V_unit**2*mp*mean_m_weight/K_B
    H = snap_H*100
    rho_bar = OmegaB*3*H**2/(8*np.pi*W.value)
    ratio_rho = snap_rho/rho_bar
    '''
    #read the position of gas
    snap_shot_gas = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','gas')
    try:
        snap_shot_star = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','star')
        snap_mass_star = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'mass','star')
    except SystemExit:
        print('no star particles now')    
    snap_SFR = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'SFR ',0)
    #the total mass distribution of snap_shot_128
    main_halo = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                         converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
    Rvir = np.array(main_halo['col12'])
    xhalo = np.array(main_halo['col6'])
    yhalo = np.array(main_halo['col7'])
    zhalo = np.array(main_halo['col8'])
    Nbins = np.array(main_halo['col37'])
    x0 = xhalo[_id_]
    y0 = yhalo[_id_]
    z0 = zhalo[_id_]
    R0 = Rvir[_id_]
    nbins0 = Nbins[_id_]
    #figure the temperature distribution of main halo 
    main_profile = pd.read_table('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%snap_z,
                                 dtype = np.float)
    r = np.array(main_profile['#r(1)'])
    r = np.abs(r)
    if _id_ == 0:
        R = r[0:nbins0]
    else:
        sum_bin = np.sum(Nbins[0:_id_])
        R = r[sum_bin : sum_bin + nbins0]
    L_R = len(R)
    
#######################section1:read out the position data of GadgetX
elif choose == 1:
    snap_z = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'redshift')
    snap_z = np.abs(snap_z)
    print('Now redshift is %.3f'%snap_z)
    snap_N = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'npartTotal')
    print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
          %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
    snap_name = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'header')
    print(snap_name)
    
    snap_mass_gas = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'mass','gas')
    snap_inE = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'U','gas') 
    snap_rho = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'RHO','gas') 
    snap_OmegaB = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'Omega0')
    print('omegaB=',snap_OmegaB)
    OmegaB = snap_OmegaB
    snap_H = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'h')
    snap_time = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'time')
    time = np.abs(snap_time)
    #PS:NE -- electron number friction is exactlly for each particle,and the same as NH.
    snap_NE = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'NE','gas')
    Ne = snap_NE
    snap_NH = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'NH','gas')
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
    snap_shot_gas = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','gas')
    try:
        snap_T = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'TEMP','gas')
    except SystemExit:
        print('no gas temperature now')
    try:
        snap_shot_star = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','star')
        snap_mass_star = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'mass','star')
    except SystemExit:
        print('no star particles now')
    snap_SFR = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'SFR ',0)
    #the total mass distribution of snap_shot_128
    main_halo = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                      converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
    Rvir = np.array(main_halo['col12'])
    xhalo = np.array(main_halo['col6'])
    yhalo = np.array(main_halo['col7'])
    zhalo = np.array(main_halo['col8'])
    Nbins = np.array(main_halo['col37'])    
    x0 = xhalo[_id_]
    y0 = yhalo[_id_]
    z0 = zhalo[_id_]
    R0 = Rvir[_id_]   
    nbins0 = Nbins[_id_]
    #figure the temperature distribution of main halo 
    main_profile = pd.read_table('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%snap_z,
                                 dtype = np.float)
    r = np.array(main_profile['#r(1)'])
    r = np.abs(r)
    if _id_ == 0:
        R = r[0:nbins0]
    else:
        sum_bin = np.sum(Nbins[0:_id_])
        R = r[sum_bin : sum_bin + nbins0]
    L_R = len(R)
else :
    print('Please selection 1(for GadgetX) or 0(for MUSIC)')
#######################section2:selection the data in goal halo
if choose == 0 or choose == 1:
    dgas = np.sqrt((snap_shot_gas[:, 0]-x0)**2 +
                   (snap_shot_gas[:, 1]-y0)**2 + (snap_shot_gas[:, 2]-z0)**2)
    ig = dgas <= R0
    inlgas = snap_shot_gas[ig, :]
    inlmass_gas = snap_mass_gas[ig]
    
    inl_rho = snap_rho[ig]
    inl_inE = snap_inE[ig]
    inl_T = snap_T[ig]
    halo_SFR = snap_SFR[ig]
    #to get the mean temperature
    mean_T_m = np.zeros(L_R,dtype = np.float)
    mean_T_n = np.zeros(L_R,dtype = np.float)
    for k in range(L_R):
        if k == 0:
            igk = dgas <= R[k]
            gasmass = snap_mass_gas[igk]
            gasT = snap_T[igk]
            mean_T_m[k] = np.sum(gasT*gasmass*10**10*Ms)/(np.sum(gasmass*10**10*Ms))
            mean_T_n[k] = np.sum(gasT)/len(gasmass)
        else:
            igk = (dgas <= R[k]) & (dgas > R[k-1])
            gasmass = snap_mass_gas[igk]
            gasT = snap_T[igk]
            mean_T_m[k] = np.sum(gasT*gasmass*10**10*Ms)/(np.sum(gasmass*10**10*Ms))
            mean_T_n[k] = np.sum(gasT)/len(gasmass)
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
        for k in range(n_r):
            if k == 0:
                igp = dgas_galaxy <= r_bcg[k]
                gasm = snap_mass_gas[igp]
                gasTg = snap_T[igp]
                T_galaxy_m[k] = np.sum(gasTg*gasm*10**10*Ms)/(np.sum(gasm*10**10*Ms))
                T_galaxy_n[k] = np.sum(gasTg)/len(gasm)
            else:
                igp = (dgas_galaxy <= r_bcg[k]) & (dgas_galaxy > r_bcg[k-1])
                gasm = snap_mass_gas[igp]
                gasTg = snap_T[igp]
                T_galaxy_m[k] = np.sum(gasTg*gasm*10**10*Ms)/(np.sum(gasm*10**10*Ms))
                T_galaxy_n[k] = np.sum(gasTg)/len(gasm)
    except ValueError:
        print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
    ##try to see the distribution of gas temperature
    #this part for MUSIC simulation
    if choose == 0:
        plt.hist2d(snap_rho,snap_SFR,bins = [100,100],normed = True,cmap = 'rainbow',
                   vmin = 1e-2,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = 'N-density')
        #handy.compare(snap_rho,snap_SFR)
        plt.xlabel(r'$\rho_{gas}[M_\odot-h^2/kpc^3]$')
        plt.ylabel(r'$SFR$')
        plt.title(r'$SFR-\rho_{gas}-MU$')
        plt.savefig('SFR-rho-gas-MU',dpi=600)
        plt.show()
        
        plt.hist2d(inl_rho,halo_SFR,bins = [100,100],normed = True,cmap = 'rainbow',
                   vmin = 1e-2,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = 'N-density')
        #handy.compare(inl_rho,halo_SFR)
        plt.xlabel(r'$\rho_{gas}[M_\odot-h^2/kpc^3]$')
        plt.ylabel(r'$SFR$')
        plt.title(r'$SFR-\rho_{gas}-mainhalo-MU$')
        plt.savefig('SFR-rho-gas-halo-MU',dpi=600)
        plt.show()
        
        plt.plot(R,mean_T_m,'r-',label = r'$\bar{T}-r-in-M$')
        plt.plot(R,mean_T_n,'b-',label = r'$\bar{T}-r-in-N$')
        plt.legend(loc = 4)
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\bar{T}-K$')
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('temperature profile mainhalo MU',dpi = 600)
        plt.show()
        
        plt.plot(r_bcg,T_galaxy_m,'r-',label = r'$\bar{T}-r-in-M$')
        plt.plot(r_bcg,T_galaxy_n,'b-',label = r'$\bar{T}-r-in-N$')        
        plt.legend(loc = 4)
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\bar{T}-K$')
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('temperatur profile BCG MU',dpi = 600)
        plt.show()
        
        plt.hist2d(np.log10(inl_rho),np.log10(inl_T),bins = [1000,1000],normed = True,
                   cmap = 'rainbow',vmin = 1e-2,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = r'$ \rho_N $')
        plt.xlabel(r'$log[\rho_g]$')
        plt.ylabel(r'$log[T]$')
        plt.title(r'$\rho_g-T$')
        plt.savefig('gas density-temperatur mainhalo MU',dpi = 600)
        plt.show()
        
        plt.hist2d(np.log10(snap_rho),np.log10(snap_T),bins = [1000,1000],normed = True,
                   cmap = 'jet',vmin = 1e-5,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = r'$ \rho_N $')
        plt.xlabel(r'$log[\rho_g]$')
        plt.ylabel(r'$log[T]$')
        plt.title(r'$\rho_g-T$')
        plt.savefig('gas density temperatur distribution MU',dpi = 600)
        plt.show() 
        
        plt.hist2d(np.log10(ratio_rho),np.log10(snap_T),bins = [1000,1000],normed = True,
                   cmap = 'jet',vmin = 1e-5,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = r'$ \rho_N $')
        plt.xlabel(r'$log[\rho_g/\rho_{bar}]$')
        plt.ylabel(r'$log[T]$')
        plt.title(r'$\rho_g-T$')
        plt.savefig('gas overdens temperatur distribution MU',dpi = 600)
        plt.show()  
    #this part for GadgetX simulation    
    if choose == 1:
        plt.hist2d(snap_rho,snap_SFR,bins = [100,100],normed = True,cmap = 'rainbow',
                   vmin = 1e-2,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = 'N-density')
        handy.compare(snap_rho,snap_SFR)
        plt.xlabel(r'$\rho_{gas}[M_\odot-h^2/kpc^3]$')
        plt.ylabel(r'$SFR$')
        plt.title(r'$SFR-\rho_{gas}-GX$')
        plt.savefig('SFR-rho-gas-GX',dpi=600)
        plt.show()
        
        plt.hist2d(inl_rho,halo_SFR,bins = [100,100],normed = True,cmap = 'rainbow',
                   vmin = 1e-2,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = 'N-density')
        handy.compare(inl_rho,halo_SFR)
        plt.xlabel(r'$\rho_{gas}[M_\odot-h^2/kpc^3]$')
        plt.ylabel(r'$SFR$')
        plt.title(r'$SFR-\rho_{gas}-mainhalo-GX$')
        plt.savefig('SFR-rho-gas-halo-GX',dpi=600)
        plt.show()
        
        plt.plot(R,mean_T_m,'r-',label = r'$\bar{T}-r-in-M$')
        plt.plot(R,mean_T_n,'b-',label = r'$\bar{T}-r-in-N$')
        plt.legend(loc = 4)
        plt.ylim(7e7,2e8)
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\bar{T}-K$')
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('temperature distribution GX',dpi = 600)
        plt.show()
        
        plt.plot(r_bcg,T_galaxy_m,'r-',label = r'$\bar{T}-r-in-M$')
        plt.plot(r_bcg,T_galaxy_n,'b-',label = r'$\bar{T}-r-in-N$')        
        plt.legend(loc = 4)
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\bar{T}-K$')
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('temperatur profile BCG GX',dpi = 600)
        plt.show()

        plt.hist2d(np.log10(inl_rho),np.log10(inl_T),bins = [1000,1000],normed = True,
                   cmap = 'rainbow',vmin = 1e-2,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = r'$ \rho_N $')
        plt.xlabel(r'$log[\rho_g]$')
        plt.ylabel(r'$log[T]$')
        plt.title(r'$\rho_g-T$')
        plt.savefig('gas density-temperatur mainhalo GX',dpi = 600)
        plt.show()
        
        plt.hist2d(np.log10(snap_rho),np.log10(snap_T),bins = [1000,1000],normed = True,
                   cmap = 'jet',vmin = 1e-5,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = r'$ \rho_N $')
        plt.xlabel(r'$log[\rho_g]$')
        plt.ylabel(r'$log[T]$')
        plt.title(r'$\rho_g-T$')
        plt.savefig('gas density temperatur distribution GX',dpi = 600)
        plt.show() 
        
        plt.hist2d(np.log10(ratio_rho),np.log10(snap_T),bins = [1000,1000],normed = True,
                   cmap = 'jet',vmin = 1e-5,vmax = 1e1,norm = mpl.colors.LogNorm())
        plt.colorbar(label = r'$ \rho_N $')
        plt.xlabel(r'$log[\rho_g/\rho_{bar}]$')
        plt.ylabel(r'$log[T]$')
        plt.title(r'$\rho_g-T$')
        plt.savefig('gas overdens temperatur distribution GX',dpi = 600)
        plt.show()                
else:
    print('Run Error!')

