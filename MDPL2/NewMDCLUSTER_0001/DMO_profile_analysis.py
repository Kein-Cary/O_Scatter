##this file using:give a redshift and the halo ID,give it's profile
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
#'D:/mask/MDPL2/NewMDCLUSTER_0001/'
##firstly,try to read the red_shift value and sort the

import glob
namestr = glob.glob('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001*AHF_halos')
#print(namestr[0])
red_shift = [float(namestr[57:-10]) for namestr in namestr]
print(red_shift)
##sort to red_shift and save as h5py file
Redshift = np.sort(red_shift)
with h5py.File('Redshift_DMO.h5','w') as f:
     f['a'] = np.array(Redshift)       
###just one time to save the redshift data
# to site the goal halo_ID and goal redshift automatically
with h5py.File('Redshift_DMO.h5') as f:
    y0 = f['a']
    Redshift = np.array(y0)
with h5py.File('main_tree_DMO.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)
##assuming now get the first one halo at z=0,and see the profile at z
goal_halo = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z0.000.AHF_halos')
goalhalo = np.array(goal_halo['#ID(1)'])
goal_id = np.int64(goalhalo[0])##set the goal halo according the data at z==0
goal_z = Redshift[0]##set the goal redshift
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

def get_data(z_value,haloid):
    z = z_value
    tr_halo = haloid#save the ID need to trace
    ###firstly,read the data at redshift,and save to review
    halo_list = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z%.3f.AHF_halos'%z)
    cNFW = np.array(halo_list['cNFW(43)'])
    Mhalo = np.array(halo_list['Mvir(4)'])
    xhalo = np.array(halo_list['Xc(6)'])
    yhalo = np.array(halo_list['Yc(7)'])
    zhalo = np.array(halo_list['Zc(8)'])
    #x,y,z: position of halo
    Rvir = np.array(halo_list['Rvir(12)'])
    #virial radius of halo
    locpeak = np.array(halo_list['r2(14)'])
    #loc_peak means the entral of potential
    Npart = np.array(halo_list['npart(5)'])
    #the number of particles of the halo
    Nbins = np.array(halo_list['nbins(37)'])
    #calculate the mass of a star particle and gasparticle 
    ARRAY = np.array([cNFW,Mhalo,xhalo,yhalo,zhalo,Rvir,locpeak,Npart,Nbins])
    #for ARRAY,each rows response one property,et. row 1--cNFW
    M = ARRAY.shape[0]
    
    #read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
    Host = np.array(halo_list['hostHalo(2)'])
    Host = np.int64(Host)
    ID = np.array(halo_list['#ID(1)'])
    ID = np.int64(ID)
    ID_save = np.array([ID, Host])
    with h5py.File('halo_ID_DMO%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
            f['a'] = np.array(ID_save)
    with h5py.File('halo_ID_DMO%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(len(ID_save)):
            f['a'][t,:] = ID_save[t,:]
            
    #get the rows of the array
    with h5py.File('halo_data_DMO%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
        f['a'] = np.array(ARRAY)
    with h5py.File('halo_data_DMO%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(M):
            f['a'][t,:] = ARRAY[t,:]
    halo_profile = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%z,dtype = np.float)
    r = np.array(halo_profile['#r(1)'])
    r = np.abs(r)#to make sure that r is positive
    n_part = np.array(halo_profile['npart(2)'])
    m_in_r = np.array(halo_profile['M_in_r(3)'])
    dens = np.array(halo_profile['dens(5)'])
    ovdens = np.array(halo_profile['ovdens(4)'])
    BRRAY = np.array([r,n_part,m_in_r,dens,ovdens])
    #for BRRAY,each rows response one property,et. row 1--r
    N = BRRAY.shape[0]
    with h5py.File('profile_data_DMO%.0f.%.3f.h5'%(tr_halo,z),'w') as f:
        f['a'] = np.array(BRRAY)
    with h5py.File('profile_data_DMO%.0f.%.3f.h5'%(tr_halo,z)) as f:
        for t in range(N):
            f['a'][t,:] = BRRAY[t,:]
    return ARRAY,BRRAY
#get_data(z_value = Redshift[0],haloid = ID[0])##test line
get_data(z_value = redshift,haloid = halo_id)
#analysis the profile 
def ana_profile(z,h_id):
    tr_halo = h_id
    r_shift = z
    #read the halo file
    with h5py.File('halo_data_DMO%.0f.%.3f.h5'%(tr_halo,r_shift)) as f:
        y1 = f['a']
        ARRAY = np.array(y1)
    '''
    cNFW = ARRAY[0,:]
    Mhalo = ARRAY[1,:]
    xhalo = ARRAY[2,:]
    yhalo = ARRAY[3,:]
    zhalo = ARRAY[4,:]
    Rvir = ARRAY[5,:]
    '''
    locpeak = ARRAY[6,:]
    Npart = ARRAY[7,:]
    Nbins = ARRAY[8,:]
    #read the profile file
    with h5py.File('profile_data_DMO%.0f.%.3f.h5'%(tr_halo,r_shift)) as f:
         y2 = f['a']
         BRRAY = np.array(y2)
    r = BRRAY[0,:]
    n_part = BRRAY[1,:]
    m_in_r = BRRAY[2,:]
    dens = BRRAY[3,:]
    ovdens = BRRAY[4,:]     
    #site the ID and then figure the profile
    with h5py.File('halo_ID_DMO%.0f.%.3f.h5'%(tr_halo,z)) as f:
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
    if _ix == 0:
        #LL = np.int(Nbins[0])#can be changed for suitting different redshift and halos ##test line
        plt.plot(r[0:LL],m_in_r[0:LL],'b-',label = r'$M_{inr}$')
        plt.legend(loc=4) 
        plt.xlabel(r'$r [kpc/h]$')
        plt.ylabel(r'$Mass [M_\odot /h]$')
        #plt.axvline(locpeak[_ix],ls='--',linewidth=0.5,color='blue')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e8,1e16)
        plt.show()
        
        #get the density
        R = r[0:LL]
        M_in = m_in_r[0:LL]
        rho_in = np.zeros(LL-1,dtype = np.float)
        for k in range(LL-1):
            dr = R[k+1]-R[k]
            rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
        #mean density and density
        plt.plot(R[0:-1],rho_in,'m-',label = r'$\rho_{tot}$')
        plt.legend(loc=1)
        plt.xlabel(r'$r [kpc/h]$')
        plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e4,1e11)
        #plt.savefig('density_profile_z_0_022.png',dpi=600)
        plt.show()
    else:
        sum_bin = np.sum(Nbins[0:_ix])
        sum_bin = np.int(sum_bin)
        #LL = np.int(Nbins[0])#can be changed for suitting different redshift and halos ##test line
        plt.plot(r[sum_bin:sum_bin+LL],m_in_r[sum_bin:sum_bin+LL],'b-',label = r'$M_{inr}$')
        plt.legend(loc=4)
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$Mass [M_\odot /h]$')
        #plt.axvline(locpeak[_ix],ls='--',linewidth=0.5,color='blue')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e8,1e16)
        plt.show()
        
        #get the density
        R = r[sum_bin:sum_bin+LL]
        M_in = m_in_r[sum_bin:sum_bin+LL]
        rho_in = np.zeros(LL-1,dtype = np.float)
        for k in range(LL-1):
            dr = R[k+1]-R[k]
            rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
        #mean density and density
        plt.plot(R[0:-1],rho_in,'m-',label = r'$\rho_{tot}$')
        plt.legend(loc=1)
        plt.xlabel(r'$r [kpc/h]$')
        plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e0,1e3)
        plt.ylim(1e4,1e11)
        #plt.savefig('density_profile_z_0_022.png',dpi=600)
        plt.show()
    return 
#ana_profile(z = redshift, h_id = ID[0])##test line
ana_profile(z = redshift, h_id = halo_id)

#next : try to figure the distrubution of the given halo
def fig_distribution(z_value,haloid):
    tr_halo = halo_id
    z = z_value
    #read the halo file
    with h5py.File('halo_data_DMO%.0f.%.3f.h5'%(tr_halo,z)) as f:
        y1 = f['a']
        ARRAY = np.array(y1)
    cNFW = ARRAY[0,:]
    Mhalo = ARRAY[1,:]
    xhalo = ARRAY[2,:]
    yhalo = ARRAY[3,:]
    zhalo =ARRAY[4,:]
    Rvir = ARRAY[5,:]
    locpeak = ARRAY[6,:]
    Npart = ARRAY[7,:]
    Nbins = ARRAY[8,:]
    #read the profile file
    with h5py.File('profile_data_DMO%.0f.%.3f.h5'%(tr_halo,z)) as f:
         y2 = f['a']
         BRRAY = np.array(y2)
    r = BRRAY[0,:]
    n_part = BRRAY[1,:]
    m_in_r = BRRAY[2,:]
    dens = BRRAY[3,:]
    ovdens = BRRAY[4,:]
    #figure
    #read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
    '''
    Host = np.array(M_halo['col2'])
    Host = np.int64(Host)
    ID = np.array(M_halo['col1'])
    #try to figure the fist one 128000000000001
    ix = Host == ID[0] ##test line
    '''
    #site the ID and then figure the profile
    with h5py.File('halo_ID_DMO%.0f.%.3f.h5'%(tr_halo,z)) as f:
        y0 = f['a'][0]
        ID = np.array(y0)
        y1 = f['a'][1]
        Host = np.array(y1)
    #site the goal ID at z!=0 from main_tree
    ip = ID == goal_progenitor
    iq = ip.tolist()
    _ix = iq.index(True)
    ix = Host == ID[_ix]
    
    #give a array to save all the position of those halos belong to 128000000000001(for others are similarly)
    position = np.zeros((len(cNFW),3),dtype = np.float)#0,1,2 coloumn respecally for : x, y, z
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
    plt.show()
    return
#fig_distribution(z_value = redshift, haloid = ID[0])##test line
fig_distribution(z_value = redshift, haloid = halo_id)
    