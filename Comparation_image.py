#this file try to see the halo image(2D)
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
idx = 0
zx = 0
# to site the goal halo_ID and goal redshift automatically
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
### get the data of GX simulation
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
    halo_profile_gx = pd.read_table('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%z,
                                    dtype = np.float)
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

def halo_image(z,hid,z_gx,hid_gx):
    tr_halo = hid
    Z = z
    #read the halo file
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_data_%.0f.%.3f.h5'%(tr_halo,Z)) as f:
        y1 = f['a']
        ARRAY = np.array(y1)
    cNFW = ARRAY[0,:]
    '''
    Mhalo = ARRAY[1,:]
    Mstar = ARRAY[2,:]
    Nstar = ARRAY[3,:]
    Mgas = ARRAY[4,:]
    Ngas = ARRAY[5,:]
    '''
    xhalo = ARRAY[6,:]
    yhalo = ARRAY[7,:]
    zhalo =ARRAY[8,:]
    Rvir = ARRAY[9,:]
    locpeak = ARRAY[10,:]
    Npart = ARRAY[11,:]
    Nbins = ARRAY[12,:]
    #read the profile file
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/profile_data_%.0f.%.3f.h5'%(tr_halo,Z)) as f:
         y2 = f['a']
         BRRAY = np.array(y2)
    '''
    r = BRRAY[0,:]
    n_part = BRRAY[1,:]
    m_in_r = BRRAY[2,:]
    m_star = BRRAY[3,:]
    m_gas = BRRAY[4,:]
    dens = BRRAY[5,:]
    ovdens = BRRAY[6,:]
    '''
    #read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
    #site the ID and then figure the profile
    with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/halo_ID_%.0f.%.3f.h5'%(tr_halo,Z)) as f:
        y3 = f['a'][0]
        ID = np.array(y3)
        y4 = f['a'][1]
        Host = np.array(y4)
    #site the goal ID at z!=0 from main_tree
    ip = ID == goal_progenitor
    iq = ip.tolist()
    _ix = iq.index(True)
    ix = Host == ID[_ix]
    #give a array to save all the position of those halos belong to 128000000000001(for others are similarly)
    position = np.zeros((len(cNFW),3),dtype = np.float)#0,1,2 coloumn respecally for : x, y, z
    
    tr_halo_gx = hid_gx
    Z_gx = z_gx
    #read the halo file
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_data_GX%.0f.%.3f.h5'%(tr_halo_gx,Z_gx)) as f:
        y1_gx = f['a']
        ARRAY_gx = np.array(y1_gx)
    cNFW_gx = ARRAY_gx[0,:]
    '''
    Mhalo_gx = ARRAY_gx[1,:]
    Mstar_gx = ARRAY_gx[2,:]
    Nstar_gx = ARRAY_gx[3,:]
    Mgas_gx = ARRAY_gx[4,:]
    Ngas_gx = ARRAY_gx[5,:]
    '''
    xhalo_gx = ARRAY_gx[6,:]
    yhalo_gx = ARRAY_gx[7,:]
    zhalo_gx = ARRAY_gx[8,:]
    Rvir_gx = ARRAY_gx[9,:]
    locpeak_gx = ARRAY_gx[10,:]
    Npart_gx = ARRAY_gx[11,:]
    Nbins_gx = ARRAY_gx[12,:]
    #read the profile file
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/profile_data_GX%.0f.%.3f.h5'%(tr_halo_gx,Z_gx)) as f:
         y2_gx = f['a']
         BRRAY_gx = np.array(y2_gx)
    '''
    r_gx = BRRAY_gx[0,:]
    n_part_gx = BRRAY_gx[1,:]
    m_in_r_gx = BRRAY_gx[2,:]
    m_star_gx = BRRAY_gx[3,:]
    m_gas_gx = BRRAY_gx[4,:]
    dens_gx = BRRAY_gx[5,:]
    ovdens_gx = BRRAY_gx[6,:]
    '''
    #read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
    #site the ID and then figure the profile
    with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/halo_ID_GX%.0f.%.3f.h5'%(tr_halo_gx,Z_gx)) as f:
        y3 = f['a'][0]
        ID_gx = np.array(y3)
        y4 = f['a'][1]
        Host_gx = np.array(y4)
    #site the goal ID at z!=0 from main_tree
    ip_gx = ID_gx == goal_progenitor_gx
    iq_gx = ip_gx.tolist()
    _ix_gx = iq_gx.index(True)
    ix_gx = Host_gx == ID_gx[_ix_gx]
    #give a array to save all the position of those halos belong to 128000000000001(for others are similarly)
    position_gx = np.zeros((len(cNFW_gx),3),dtype = np.float)#0,1,2 coloumn respecally for : x, y, z
    plt.figure()
    plt.subplot(1,2,1)
    for k in range(len(ix)):
        if k == _ix :
            position[k,:] = np.array([xhalo[k],yhalo[k],zhalo[k]])
        elif ix[k] == True :
            position[k,:] = np.array([xhalo[k],yhalo[k],zhalo[k]])
        else:
            position[k,:] = np.array([np.inf,np.inf,np.inf])
    for k in range(len(ID)):
        iy = np.inf in position[k,:]
        if iy == False :
           #hsc.circles(xhalo[k],yhalo[k],s=Rvir[k],c = 'g', alpha = 0.2)
           hsc.circles(xhalo[k],zhalo[k],s=Rvir[k],c = 'g', alpha = 0.2)
           #hsc.circles(yhalo[k],zhalo[k],s=Rvir[k],c = 'g', alpha = 0.2)
    plt.title(r'$MUSIC$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    #plt.show()
    plt.subplot(1,2,2)
    for k in range(len(ix_gx)):
        if k == _ix_gx :
            position_gx[k,:] = np.array([xhalo_gx[k],yhalo_gx[k],zhalo_gx[k]])
        elif ix_gx[k] == True :
            position_gx[k,:] = np.array([xhalo_gx[k],yhalo_gx[k],zhalo_gx[k]])
        else:
            position_gx[k,:] = np.array([np.inf,np.inf,np.inf])
    for k in range(len(ID_gx)):
        iy_gx = np.inf in position_gx[k,:]
        if iy_gx == False :
           #hsc.circles(xhalo_gx[k],yhalo_gx[k],s=Rvir_gx[k],c = 'g', alpha = 0.2)
           hsc.circles(xhalo_gx[k],zhalo_gx[k],s=Rvir_gx[k],c = 'g', alpha = 0.2)
           #hsc.circles(yhalo_gx[k],zhalo_gx[k],s=Rvir_gx[k],c = 'g', alpha = 0.2)
    plt.title(r'$GX$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    #plt.show()
    plt.tight_layout()
    plt.savefig('Image halo',dpi=600)
    plt.show()
    return
halo_image(z=redshift,hid=halo_id,z_gx=redshift_gx,hid_gx=halo_id_gx)