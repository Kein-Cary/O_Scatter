##This fiel for reading data of AHF
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
#'D:/mask/MDPL2/NewMDCLUSTER_0001/'
####part 0:read the redshift,to show the result in different z
with h5py.File('Redshift_DMO.h5') as f:
    y0 = f['a']
    z = np.array(y0)
with h5py.File('main_tree_DMO.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)
##assuming now get the first one halo at z=0,and see the profile at z
goal_halo = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z0.000.AHF_halos')
goalhalo = np.array(goal_halo['#ID(1)'])
goal_id = np.int64(goalhalo[0])##set the goal halo
goal_z = z[0]##set the goal redshift
IX = z == goal_z
IY = IX.tolist()
No_site_z = IY.index(True)
ia = main_tree[:,0] == goal_id
ib = ia.tolist()
No_site_id = ib.index(True)
#get the goal main progenitor
goal_progenitor = main_tree[No_site_id,No_site_z]

###############part1
M_halo = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z%.3f.AHF_halos'%z[No_site_z])
# try to get the base property of halo 128000000000001
cNFW = np.array(M_halo['cNFW(43)'])
Mhalo = np.array(M_halo['Mvir(4)'])
xhalo = np.array(M_halo['Xc(6)'])
yhalo = np.array(M_halo['Yc(7)'])
zhalo = np.array(M_halo['Zc(8)'])
#x,y,z: position of halo
Rvir = np.array(M_halo['Rvir(12)'])
#virial radius of halo
locpeak = np.array(M_halo['r2(14)'])
#loc_peak means the entral of potential
Npart = np.array(M_halo['npart(5)'])
#the number of particles of the halo
Nbins = np.array(M_halo['nbins(37)'])
#read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
Host = np.array(M_halo['hostHalo(2)'])
Host = np.int64(Host)
ID = np.array(M_halo['#ID(1)'])
ID = np.int64(ID)

#site the goal ID at z!=0 from main_tree
ip = ID == goal_progenitor
iq = ip.tolist()
_ix = iq.index(True)
ix = Host == ID[_ix]

#try to figure the fist one 128000000000001
#ix = Host == ID[0]
#give a array to save all the position of those halos belong to 128000000000001(for others are similarly)
position = np.zeros((len(cNFW),3),dtype = np.float)#0,1,2 coloumn respecally for : x, y, z
for k in range(len(ix)):
    if k == _ix:
       position[k,:] = np.array([xhalo[_ix],yhalo[_ix],zhalo[_ix]]) 
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
################part2
h_profile = pd.read_table('D:/mask/MDPL2/NewMDCLUSTER_0001/MDPL2-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%z[No_site_z],dtype = np.float)
r = np.array(h_profile['#r(1)'])
r = np.abs(r)#to make sure that r is positive
n_part = np.array(h_profile['npart(2)'])
m_in_r = np.array(h_profile['M_in_r(3)'])
dens = np.array(h_profile['dens(5)'])
ovdens = np.array(h_profile['ovdens(4)'])
###############part3
##next : show the properties of first halo
#N = np.int(Nbins[0])#can be changed for suitting different redshift and halos

N = np.int(Nbins[_ix])
if _ix==0 :
    plt.plot(r[0:N],m_in_r[0:N],'b-',label = r'$M_{inr}$')
    plt.legend(loc=4)
    plt.xlabel(r'$r [kpc/h]$')
    plt.ylabel(r'$Mass [M_\odot /h]$')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1e0,1e3)
    plt.ylim(1e8,1e16)
    plt.show()
    
    #get the density
    R = r[0:N]
    M_in = m_in_r[0:N]
    rho_in = np.zeros(N-1,dtype = np.float)
    for k in range(N-1):
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
    plt.show()
else:
    sum_bin = np.sum(Nbins[0:_ix])
    plt.plot(r[sum_bin:sum_bin+N],m_in_r[sum_bin:sum_bin+N],'b-',label = r'$M_{inr}$')    
    plt.legend(loc=4)
    plt.xlabel(r'$r [kpc/h]$')
    plt.ylabel(r'$Mass [M_\odot /h]$')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1e0,1e3)
    plt.ylim(1e8,1e16)
    plt.show()
    
    #get the density
    R = r[sum_bin:sum_bin+N]
    M_in = m_in_r[sum_bin:sum_bin+N]
    rho_in = np.zeros(N-1,dtype = np.float)
    for k in range(N-1):
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
    #plt.savefig('density_profile_DMO',dpi=600)
    plt.show()
#try to count the particles of 



