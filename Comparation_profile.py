#this file try to comparation the mass profile and density profile between GX and MUSIC simulation
#for give halo at z==0 and given redshift z~=0
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
#import glob
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
#set the goal halo and goal redshift
idx = 0
zx = 0
#to site the goal halo_ID and goal redshift automatically
with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/Redshift.h5') as f:
    y0 = f['a']
    Redshift = np.array(y0)
with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/main_tree.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)
##assuming now get the first one halo at z=0,and see the profile at z
goal_halo = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
goalhalo = np.array(goal_halo['col1'])
goal_id = np.int64(goalhalo[idx])##set the goal halo according the data at z==0
goal_z = Redshift[zx]##set the goal redshift
IX = Redshift == goal_z
IY = IX.tolist()
No_site_z = IY.index(True)
ia = main_tree[:,0] == goal_id
ib = ia.tolist()
No_site_id = ib.index(True)
#get the goal main progenitor
goal_progenitor = main_tree[No_site_id,No_site_z]
# use the site number to find the result
redshift = goal_z
halo_id = goal_id

#get data for GX simulation
with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/Redshift_GX.h5') as f:
    y2 = f['a']
    Redshift_gx = np.array(y2)
with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/main_tree_GX.h5') as f:
    y3 = f['a']
    main_tree_gx = np.array(y3)
##assuming now get the first one halo at z=0,and see the profile at z
goal_halo_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
goalhalo_gx = np.array(goal_halo_gx['col1'])
goal_id_gx = np.int64(goalhalo_gx[idx])##set the goal halo according the data at z==0
goal_z_gx = Redshift_gx[zx]##set the goal redshift
IX_gx = Redshift_gx == goal_z_gx
IY_gx = IX_gx.tolist()
No_site_z_gx = IY_gx.index(True)
ia_gx = main_tree_gx[:,0] == goal_id_gx
ib_gx = ia_gx.tolist()
No_site_id_gx = ib_gx.index(True)
#get the goal main progenitor
goal_progenitor_gx = main_tree_gx[No_site_id_gx,No_site_z_gx]
# use the site number to find the result
redshift_gx = goal_z_gx
halo_id_gx = goal_id_gx

def get_data(z_value,haloid):
    z = z_value
    tr_halo = haloid#save the ID need to trace
    ###firstly,read the data at redshift,and save to review
    halo_list = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%z,
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
    cNFW = np.array(halo_list['col43'])
    Mhalo = np.array(halo_list['col4'])
    Mstar = np.array(halo_list['col65'])
    Nstar = np.array(halo_list['col64'])
    #star mass and number
    Mgas = np.array(halo_list['col45'])
    Ngas = np.array(halo_list['col44'])
    #gas mass and number
    xhalo = np.array(halo_list['col6'])
    yhalo = np.array(halo_list['col7'])
    zhalo = np.array(halo_list['col8'])
    #x,y,z: position of halo
    Rvir = np.array(halo_list['col12'])
    #virial radius of halo
    locpeak = np.array(halo_list['col14'])
    #loc_peak means the entral of potential
    Npart = np.array(halo_list['col5'])
    #the number of particles of the halo
    Nbins = np.array(halo_list['col37'])
    #calculate the mass of a star particle and gasparticle 
    ARRAY = np.array([cNFW,Mhalo,Mstar,Nstar,Mgas,Ngas,xhalo,yhalo,zhalo,Rvir,locpeak,Npart,Nbins])
    #for ARRAY,each rows response one property,et. row 1--cNFW
    M = ARRAY.shape[0]
    #read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
    Host = np.array(halo_list['col2'])
    Host = np.int64(Host)
    ID = np.array(halo_list['col1'])
    ID = np.int64(ID)
    ID_save = np.array([ID, Host])
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_ID_%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
            f['a'] = np.array(ID_save)
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_ID_%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(len(ID_save)):
            f['a'][t,:] = ID_save[t,:]
            
    #get the rows of the array
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_data_%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
        f['a'] = np.array(ARRAY)
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_data_%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(M):
            f['a'][t,:] = ARRAY[t,:]
    halo_profile = pd.read_table('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%z,dtype = np.float)
    r = np.array(halo_profile['#r(1)'])
    r = np.abs(r)#to make sure that r is positive
    n_part = np.array(halo_profile['npart(2)'])
    m_in_r = np.array(halo_profile['M_in_r(3)'])
    m_star = np.array(halo_profile['M_star(26)'])
    m_gas = np.array(halo_profile['M_gas(25)'])
    dens = np.array(halo_profile['dens(5)'])
    ovdens = np.array(halo_profile['ovdens(4)'])
    BRRAY = np.array([r,n_part,m_in_r,m_star,m_gas,dens,ovdens])
    #for BRRAY,each rows response one property,et. row 1--r
    N = BRRAY.shape[0]
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/profile_data_%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
        f['a'] = np.array(BRRAY)
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/profile_data_%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(N):
            f['a'][t,:] = BRRAY[t,:]       
    return ARRAY,BRRAY
get_data(z_value = redshift,haloid = halo_id)

#get data of GX simulation
def get_data_gx(z_value,haloid):
    z = z_value
    tr_halo = haloid#save the ID need to trace
    ###firstly,read the data at redshift,and save to review
    halo_list_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%z,
                converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
    cNFW_gx = np.array(halo_list_gx['col43'])
    Mhalo_gx = np.array(halo_list_gx['col4'])
    Mstar_gx = np.array(halo_list_gx['col65'])
    Nstar_gx = np.array(halo_list_gx['col64'])
    #star mass and number
    Mgas_gx = np.array(halo_list_gx['col45'])
    Ngas_gx = np.array(halo_list_gx['col44'])
    #gas mass and number
    xhalo_gx = np.array(halo_list_gx['col6'])
    yhalo_gx = np.array(halo_list_gx['col7'])
    zhalo_gx = np.array(halo_list_gx['col8'])
    #x,y,z: position of halo
    Rvir_gx = np.array(halo_list_gx['col12'])
    #virial radius of halo
    locpeak_gx = np.array(halo_list_gx['col14'])
    #loc_peak means the entral of potential
    Npart_gx = np.array(halo_list_gx['col5'])
    #the number of particles of the halo
    Nbins_gx = np.array(halo_list_gx['col37'])
    #calculate the mass of a star particle and gasparticle 
    ARRAY_gx = np.array([cNFW_gx,Mhalo_gx,Mstar_gx,Nstar_gx,Mgas_gx,Ngas_gx,
                         xhalo_gx,yhalo_gx,zhalo_gx,Rvir_gx,locpeak_gx,Npart_gx,Nbins_gx])
    #for ARRAY,each rows response one property,et. row 1--cNFW
    M_gx = ARRAY_gx.shape[0]
    
    #read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
    Host_gx = np.array(halo_list_gx['col2'])
    Host_gx = np.int64(Host_gx)
    ID_gx = np.array(halo_list_gx['col1'])
    ID_gx = np.int64(ID_gx)
    ID_save_gx = np.array([ID_gx, Host_gx])
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_ID_GX%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
            f['a'] = np.array(ID_save_gx)
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_ID_GX%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(len(ID_save_gx)):
            f['a'][t,:] = ID_save_gx[t,:]
            
    #get the rows of the array
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_data_GX%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
        f['a'] = np.array(ARRAY_gx)
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_data_GX%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(M_gx):
            f['a'][t,:] = ARRAY_gx[t,:]
    halo_profile_gx = pd.read_table('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%z,dtype = np.float)
    r_gx = np.array(halo_profile_gx['#r(1)'])
    r_gx = np.abs(r_gx)#to make sure that r is positive
    n_part_gx = np.array(halo_profile_gx['npart(2)'])
    m_in_r_gx = np.array(halo_profile_gx['M_in_r(3)'])
    m_star_gx = np.array(halo_profile_gx['M_star(26)'])
    m_gas_gx = np.array(halo_profile_gx['M_gas(25)'])
    dens_gx = np.array(halo_profile_gx['dens(5)'])
    ovdens_gx = np.array(halo_profile_gx['ovdens(4)'])
    BRRAY_gx = np.array([r_gx,n_part_gx,m_in_r_gx,m_star_gx,m_gas_gx,dens_gx,ovdens_gx])
    #for BRRAY,each rows response one property,et. row 1--r
    N_gx = BRRAY_gx.shape[0]
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/profile_data_GX%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
        f['a'] = np.array(BRRAY_gx)
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/profile_data_GX%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(N_gx):
            f['a'][t,:] = BRRAY_gx[t,:]
    return ARRAY_gx,BRRAY_gx
get_data_gx(z_value = redshift_gx,haloid = halo_id_gx)

#Analysis the profile 
def ana_profile(z,h_id,z_gx,h_id_gx):
    tr_halo = h_id
    r_shift = z
    #read the halo file
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_data_%.0f.%.3f.h5'%(tr_halo,r_shift)) as f:
        y_1 = f['a']
        ARRAY = np.array(y_1)
    locpeak = ARRAY[10,:]
    Npart = ARRAY[11,:]
    Nbins = ARRAY[12,:]
    #read the profile file
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/profile_data_%.0f.%.3f.h5'%(tr_halo,r_shift)) as f:
         y_2 = f['a']
         BRRAY = np.array(y_2)
    r = BRRAY[0,:]
    n_part = BRRAY[1,:]
    m_in_r = BRRAY[2,:]
    m_star = BRRAY[3,:]
    m_gas = BRRAY[4,:]
    dens = BRRAY[5,:]
    ovdens = BRRAY[6,:]
    #site the ID and then figure the profile
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_ID_%.0f.%.3f.h5'%(tr_halo,r_shift)) as f:
        y3 = f['a'][0]
        ID = np.array(y3)
        y4 = f['a'][1]
        Host = np.array(y4)
    #site the goal ID at z!=0 from main_tree
    ip = ID == goal_progenitor
    iq = ip.tolist()
    _ix = iq.index(True)
    ix = Host == ID[_ix]
    ##next : show the properties of goal halo
    LL = np.int(Nbins[_ix])
    
    #get the data of GX simulation
    tr_halo_gx = h_id_gx
    r_shift_gx = z_gx
    #read the halo file
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_data_GX%.0f.%.3f.h5'%(tr_halo_gx,r_shift_gx)) as f:
        y_1_gx = f['a']
        ARRAY_gx = np.array(y_1_gx)
    locpeak_gx = ARRAY_gx[10,:]
    Npart_gx = ARRAY_gx[11,:]
    Nbins_gx = ARRAY_gx[12,:]
    #read the profile file
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/profile_data_GX%.0f.%.3f.h5'%(tr_halo_gx,r_shift_gx)) as f:
         y_2_gx = f['a']
         BRRAY_gx = np.array(y_2_gx)
    r_gx = BRRAY_gx[0,:]
    n_part_gx = BRRAY_gx[1,:]
    m_in_r_gx = BRRAY_gx[2,:]
    m_star_gx = BRRAY_gx[3,:]
    m_gas_gx = BRRAY_gx[4,:]
    dens_gx = BRRAY_gx[5,:]
    ovdens_gx = BRRAY_gx[6,:]
    #site the ID and then figure the profile
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_ID_GX%.0f.%.3f.h5'%(tr_halo_gx,r_shift_gx)) as f:
        y_3_gx = f['a'][0]
        ID_gx = np.array(y_3_gx)
        y_4_gx = f['a'][1]
        Host_gx = np.array(y_4_gx)
    #site the goal ID at z!=0 from main_tree
    ip_gx = ID_gx == goal_progenitor_gx
    iq_gx = ip_gx.tolist()
    _ix_gx = iq_gx.index(True)
    ix_gx = Host_gx == ID_gx[_ix_gx]
    ##next : show the properties of goal halo
    LL_gx = np.int(Nbins_gx[_ix_gx])
    #figure the profile
    plt.figure()
    if _ix == 0:
        m_b = m_gas[0:LL] + m_star[0:LL]
        m_h = m_in_r[0:LL] - m_b
        plt.plot(r[0:LL],m_in_r[0:LL],'b-',label = r'$M_{inr}$')
        plt.plot(r[0:LL],m_star[0:LL],'r-',label = r'$M_\ast$')
        plt.plot(r[0:LL],m_gas[0:LL],'g-',label = r'$M_{gas}$')
        plt.plot(r[0:LL],m_b,'m-',label = r'$M_b$')
        plt.plot(r[0:LL],m_h,'k-',label = r'$M_d$')
        #plt.legend(loc=4) 
        #plt.xlabel(r'$r [kpc/h]$')
        #plt.ylabel(r'$Mass [M_\odot /h]$')
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.xlim(1e0,1e3)
        #plt.ylim(1e8,1e16)
        #plt.show()
    else:
        sum_bin = np.sum(Nbins[0:_ix])
        sum_bin = np.int(sum_bin)
        m_b = m_gas[sum_bin:sum_bin+LL] + m_star[sum_bin:sum_bin+LL]
        m_h = m_in_r[sum_bin:sum_bin+LL] - m_b
        plt.plot(r[sum_bin:sum_bin+LL],m_in_r[sum_bin:sum_bin+LL],'b-',label = r'$M_{inr}$')
        plt.plot(r[sum_bin:sum_bin+LL],m_star[sum_bin:sum_bin+LL],'r-',label = r'$M_\ast$')
        plt.plot(r[sum_bin:sum_bin+LL],m_gas[sum_bin:sum_bin+LL],'g-',label = r'$M_{gas}$')
        plt.plot(r[sum_bin:sum_bin+LL],m_b,'m-',label = r'$M_b$')
        plt.plot(r[sum_bin:sum_bin+LL],m_h,'k-',label = r'$M_d$')
        #plt.legend(loc=4)
        #plt.xlabel(r'$r-kpc/h$')
        #plt.ylabel(r'$Mass [M_\odot /h]$')
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.xlim(1e0,1e3)
        #plt.ylim(1e8,1e16)
        #plt.show() 
    if _ix_gx == 0:
        m_b_gx = m_gas_gx[0:LL_gx] + m_star_gx[0:LL_gx]
        m_h_gx = m_in_r_gx[0:LL_gx] - m_b_gx
        plt.plot(r_gx[0:LL_gx],m_in_r_gx[0:LL_gx],'b--',label = r'$M_{inr} - GX$')
        plt.plot(r_gx[0:LL_gx],m_star_gx[0:LL_gx],'r--',label = r'$M_\ast -GX$')
        plt.plot(r_gx[0:LL_gx],m_gas_gx[0:LL_gx],'g--',label = r'$M_{gas} -GX$')
        plt.plot(r_gx[0:LL_gx],m_b_gx,'m--',label = r'$M_b -GX$')
        plt.plot(r_gx[0:LL_gx],m_h_gx,'k--',label = r'$M_d -GX$')
        plt.legend(bbox_to_anchor=(0., 1., 1.15, 0.), loc=2,
           ncol=2, fontsize=7.5, borderaxespad=0.) 
        plt.xlabel(r'$r [kpc/h]$')
        plt.ylabel(r'$Mass [M_\odot /h]$')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e8,1e16)
        plt.title('Mass profile')
        plt.savefig('Mass profile',dpi=600)
        plt.show()
    else:
        sum_bin_gx = np.sum(Nbins_gx[0:_ix_gx])
        sum_bin_gx = np.int(sum_bin_gx)
        m_b_gx = m_gas_gx[sum_bin_gx:sum_bin_gx+LL_gx] + m_star_gx[sum_bin_gx:sum_bin_gx+LL_gx]
        m_h_gx = m_in_r_gx[sum_bin_gx:sum_bin_gx+LL_gx] - m_b_gx
        plt.plot(r_gx[sum_bin_gx:sum_bin_gx+LL_gx], m_in_r_gx[sum_bin_gx:sum_bin_gx+LL_gx],'b--',
                 label = r'$M_{inr} -GX$')
        plt.plot(r_gx[sum_bin_gx:sum_bin_gx+LL_gx], m_star_gx[sum_bin_gx:sum_bin_gx+LL_gx],'r--',
                 label = r'$M_\ast -GX$')
        plt.plot(r_gx[sum_bin_gx:sum_bin_gx+LL_gx], m_gas_gx[sum_bin_gx:sum_bin_gx+LL_gx],'g--',
                 label = r'$M_{gas} -GX$')
        plt.plot(r_gx[sum_bin_gx:sum_bin_gx+LL_gx],m_b_gx,'m--',label = r'$M_b -GX$')
        plt.plot(r_gx[sum_bin_gx:sum_bin_gx+LL_gx],m_h_gx,'k--',label = r'$M_d -GX$')
        plt.legend(bbox_to_anchor=(0., 1., 1.15, 0.), loc=2,
           ncol=2, fontsize=7.5, borderaxespad= 0.) 
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$Mass [M_\odot /h]$')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e8,1e16)
        plt.title('Mass profile')
        plt.savefig('Mass_profile',dpi=600)
        plt.show()
    plt.figure()
    if _ix == 0:
        #get the density
        R = r[0:LL]
        M_in = m_in_r[0:LL]
        M_in_gas = m_gas[0:LL]
        M_star = m_star[0:LL]
        rho_in = np.zeros(LL-1,dtype = np.float)
        rho_h = np.zeros(LL-1,dtype = np.float)
        rho_g = np.zeros(LL-1,dtype = np.float)
        rho_s = np.zeros(LL-1,dtype = np.float)
        for k in range(LL-1):
            dr = R[k+1]-R[k]
            rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
            rho_h[k] = (m_h[k+1]-m_h[k])/(4*np.pi*R[k]**2*dr)
            rho_g[k] = (M_in_gas[k+1]-M_in_gas[k])/(4*np.pi*R[k]**2*dr)
            rho_s[k] = (M_star[k+1]-M_star[k])/(4*np.pi*R[k]**2*dr)
        #mean density and density
        plt.plot(R[0:-1],rho_in,'m-',label = r'$\rho_{tot}$')
        plt.plot(R[0:-1],rho_h,'k-',label = r'$\rho_d$')
        plt.plot(R[0:-1],rho_g,'g-',label = r'$\rho_{gas}$')
        plt.plot(R[0:-1],rho_s,'r-',label = r'$\rho_{star}$')
        #plt.legend(loc=1)
        #plt.xlabel(r'$r [kpc/h]$')
        #plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.xlim(1e0,1e3)
        #plt.ylim(1e4,1e11)
        #plt.savefig('density_profile',dpi=600)
        #plt.show()
    else:
        #get the density
        R = r[sum_bin:sum_bin+LL]
        M_in = m_in_r[sum_bin:sum_bin+LL]
        M_in_gas = m_gas[sum_bin:sum_bin+LL]
        M_star = m_star[sum_bin:sum_bin+LL]
        rho_in = np.zeros(LL-1,dtype = np.float)
        rho_h = np.zeros(LL-1,dtype = np.float)
        rho_g = np.zeros(LL-1,dtype = np.float)
        rho_s = np.zeros(LL-1,dtype = np.float)
        for k in range(LL-1):
            dr = R[k+1]-R[k]
            rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
            rho_h[k] = (m_h[k+1]-m_h[k])/(4*np.pi*R[k]**2*dr)
            rho_g[k] = (M_in_gas[k+1]-M_in_gas[k])/(4*np.pi*R[k]**2*dr)
            rho_s[k] = (M_star[k+1]-M_star[k])/(4*np.pi*R[k]**2*dr)
        plt.plot(R[0:-1],rho_in,'m-',label = r'$\rho_{tot}$')
        plt.plot(R[0:-1],rho_h,'k-',label = r'$\rho_d$')
        plt.plot(R[0:-1],rho_g,'g-',label = r'$\rho_{gas}$')
        plt.plot(R[0:-1],rho_s,'r-',label = r'$\rho_{star}$')
        #plt.legend(loc=1)
        #plt.xlabel(r'$r [kpc/h]$')
        #plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.xlim(1e0,1e3)
        #plt.ylim(1e4,1e11)
        #plt.savefig('density_profile',dpi=600)
        #plt.show()
    if _ix_gx == 0:
        R_gx = r_gx[0:LL_gx]
        M_in_gx = m_in_r_gx[0:LL_gx]
        M_in_gas_gx = m_gas_gx[0:LL_gx]
        M_star_gx = m_star_gx[0:LL_gx]
        rho_in_gx = np.zeros(LL_gx-1,dtype = np.float)
        rho_h_gx = np.zeros(LL_gx-1,dtype = np.float)
        rho_g_gx = np.zeros(LL_gx-1,dtype = np.float)
        rho_s_gx = np.zeros(LL_gx-1,dtype = np.float)
        for k in range(LL_gx-1):
            dr_gx = R_gx[k+1]-R_gx[k]
            rho_in_gx[k] = (M_in_gx[k+1]-M_in_gx[k])/(4*np.pi*R_gx[k]**2*dr_gx)
            rho_h_gx[k] = (m_h_gx[k+1]-m_h_gx[k])/(4*np.pi*R_gx[k]**2*dr_gx)
            rho_g_gx[k] = (M_in_gas_gx[k+1]-M_in_gas_gx[k])/(4*np.pi*R_gx[k]**2*dr_gx)
            rho_s_gx[k] = (M_star_gx[k+1]-M_star_gx[k])/(4*np.pi*R_gx[k]**2*dr)
        #mean density and density
        plt.plot(R_gx[0:-1],rho_in_gx,'m--',label = r'$\rho_{tot} -GX$')
        plt.plot(R_gx[0:-1],rho_h_gx,'k--',label = r'$\rho_d -GX$')
        plt.plot(R_gx[0:-1],rho_g_gx,'g--',label = r'$\rho_{gas} -GX$')
        plt.plot(R_gx[0:-1],rho_s_gx,'r--',label = r'$\rho_{star} -GX$')
        plt.legend(loc=1)
        plt.xlabel(r'$r [kpc/h]$')
        plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e4,1e11)
        plt.title('density profile')
        plt.savefig('density_profile',dpi=600)
        plt.show() 
    else:
        R_gx = r_gx[sum_bin_gx:sum_bin_gx+LL_gx]
        M_in_gx = m_in_r_gx[sum_bin_gx:sum_bin_gx+LL_gx]
        M_in_gas_gx = m_gas_gx[sum_bin_gx:sum_bin_gx+LL_gx]
        M_star_gx = m_star_gx[sum_bin_gx:sum_bin_gx+LL_gx]
        rho_in_gx = np.zeros(LL_gx-1,dtype = np.float)
        rho_h_gx = np.zeros(LL_gx-1,dtype = np.float)
        rho_g_gx = np.zeros(LL_gx-1,dtype = np.float)
        rho_s_gx = np.zeros(LL_gx-1,dtype = np.float)
        for k in range(LL_gx-1):
            dr_gx = R_gx[k+1]-R_gx[k]
            rho_in_gx[k] = (M_in_gx[k+1]-M_in_gx[k])/(4*np.pi*R_gx[k]**2*dr_gx)
            rho_h_gx[k] = (m_h_gx[k+1]-m_h_gx[k])/(4*np.pi*R_gx[k]**2*dr_gx)
            rho_g_gx[k] = (M_in_gas_gx[k+1]-M_in_gas_gx[k])/(4*np.pi*R_gx[k]**2*dr_gx)
            rho_s_gx[k] = (M_star_gx[k+1]-M_star_gx[k])/(4*np.pi*R_gx[k]**2*dr)
        plt.plot(R_gx[0:-1],rho_in_gx,'m--',label = r'$\rho_{tot} -GX$')
        plt.plot(R_gx[0:-1],rho_h_gx,'k--',label = r'$\rho_d -GX$')
        plt.plot(R_gx[0:-1],rho_g_gx,'g--',label = r'$\rho_{gas} -GX$')
        plt.plot(R_gx[0:-1],rho_s_gx,'r--',label = r'$\rho_{star} -GX$')
        plt.legend(loc=1)
        plt.xlabel(r'$r [kpc/h]$')
        plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e4,1e11)
        plt.title('density profile')
        plt.savefig('density_profile',dpi=600)
        plt.show()
    return 
ana_profile(z = redshift, h_id = halo_id, z_gx = redshift_gx, h_id_gx = halo_id_gx)