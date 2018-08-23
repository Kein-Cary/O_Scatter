##This fiel for reading data of AHF
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
sys.path.insert(0,'D:/mask')
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
###############part1
M_halo = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
                  ## this line try to set the datatype of the first coloum
# try to get the base property of halo 128000000000001
#test line
#M_halo = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z1.321.AHF_halos',
#                 converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
cNFW = np.array(M_halo['col43'])
Mhalo = np.array(M_halo['col4'])
Mstar = np.array(M_halo['col65'])
Nstar = np.array(M_halo['col64'])
#star mass and number
Mgas = np.array(M_halo['col45'])
Ngas = np.array(M_halo['col44'])
#gas mass and number
xhalo = np.array(M_halo['col6'])
yhalo = np.array(M_halo['col7'])
zhalo = np.array(M_halo['col8'])
#x,y,z: position of halo
Rvir = np.array(M_halo['col12'])
#virial radius of halo
locpeak = np.array(M_halo['col14'])
#loc_peak means the entral of potential
Npart = np.array(M_halo['col5'])
#the number of particles of the halo
Nbins = np.array(M_halo['col37'])
#read the hosthalo value,and the ID value,and then compare to find the number of halo belong to a give ID
Host = np.array(M_halo['col2'])
Host = np.int64(Host)
ID = np.array(M_halo['col1'])
#try to figure the fist one 128000000000001
_id_ = 0
ix = Host == ID[_id_]
#give a array to save all the position of those halos belong to 128000000000001(for others are similarly)
position = np.zeros((len(cNFW),3),dtype = np.float)#0,1,2 coloumn respecally for : x, y, z
#save the link Rvir for those halo
for k in range(len(ix)):
    if k == _id_:
       position[k,:] = np.array([xhalo[_id_],yhalo[_id_],zhalo[_id_]]) 
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
#plt.savefig('halo_view_y_z_0_022.png',dpi=600)
plt.show()
################part2
h_profile = pd.read_table('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_profiles',dtype = np.float)
#test line
#h_profile = pd.read_table('GadgetMUSIC-NewMDCLUSTER_0001.z1.321.AHF_profiles',dtype = np.float)
r = np.array(h_profile['#r(1)'])
r = np.abs(r)#to make sure that r is positive
n_part = np.array(h_profile['npart(2)'])
m_in_r = np.array(h_profile['M_in_r(3)'])
m_star = np.array(h_profile['M_star(26)'])
m_gas = np.array(h_profile['M_gas(25)'])
dens = np.array(h_profile['dens(5)'])
ovdens = np.array(h_profile['ovdens(4)'])
dm = m_in_r - m_star - m_gas
###############part3
#color for particles,0:gas--green,1:DM--M,4:star--r;the total mass:b
color = ['g','M','r','b']
##next : show the properties of first halo
N = np.int(Nbins[_id_])#can be changed for suitting different redshift and halos
if _id_ == 0:
    plt.plot(r[0:N],m_in_r[0:N],'b-',label = r'$M_{inr}$')
    plt.plot(r[0:N],m_star[0:N],'r-',label = r'$M_\ast$')
    plt.plot(r[0:N],m_gas[0:N],'g-',label = r'$M_{gas}$')
    
    m_b = m_gas[0:N] + m_star[0:N]
    plt.plot(r[0:N],m_b,'m-',label = r'$M_b$')
    m_h = m_in_r[0:N] - m_b
    plt.plot(r[0:N],m_h,'k-',label = r'$M_d$')
    
    plt.legend(loc=4)
    plt.xlabel(r'$r [kpc/h]$')
    plt.ylabel(r'$Mass [M_\odot /h]$')
    #plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1e0,1e3)
    plt.ylim(1e8,1e16)
    #plt.savefig('mass_profile_z_0_022.png',dpi=600)
    plt.show()
    
    #get the density
    R = r[0:N]
    M_in = m_in_r[0:N]
    M_in_gas = m_gas[0:N]
    rho_in = np.zeros(N-1,dtype = np.float)
    rho_h = np.zeros(N-1,dtype = np.float)
    rho_g = np.zeros(N-1,dtype = np.float)
    for k in range(N-1):
        dr = R[k+1]-R[k]
        rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
        rho_h[k] = (m_h[k+1]-m_h[k])/(4*np.pi*R[k]**2*dr)
        rho_g[k] = (M_in_gas[k+1]-M_in_gas[k])/(4*np.pi*R[k]**2*dr)
    #mean density and density
    plt.plot(R[0:-1],rho_in,'m-',label = r'$\rho_{tot}$')
    plt.plot(R[0:-1],rho_h,'k-',label = r'$\rho_d$')
    plt.plot(R[0:-1],rho_g,'g-',label = r'$\rho_{gas}$')
    plt.legend(loc=1)
    plt.xlabel(r'$r [kpc/h]$')
    plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1e0,1e3)
    plt.ylim(1e4,1e11)
    #plt.savefig('density_profile_z_0_022.png',dpi=600)
    plt.show()
    '''
    plt.plot(r[0:N],dens[0:N],'g--',label = r'$\rho$')
    plt.plot(r[0:N],ovdens[0:N],'r--',label = r'$\Delta \rho$')
    plt.legend(loc=4)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$M_\odot h^2 {kpc^3}$')
    #plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    #plt.savefig('density_profile_z_0_022.png',dpi=600)
    plt.show()
    #nbins:use this to identical the particles' rows belong to the halo in the *AHF_particles file
    m_p_star = Mstar[0]/Nstar[0]
    m_p_gas = Mgas[0]/Ngas[0]
    plt.plot(r[0:N],n_part[0:N],'g',label = 'Num-particle')
    plt.legend(loc = 4)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'N')
    #plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    #plt.savefig('Number_z_0_022.png',dpi=600)
    plt.show()
    
    plt.plot(r[0:N],n_part[0:N]/(4*np.pi*r[0:N]**3/3),'r',label = 'Num-density')
    plt.legend(loc = 4)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$N h^3/kpc^3$')
    #plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    #plt.savefig('Number_density_z_0_022')
    plt.show()
    '''
else:
    sum_bin = np.sum(Nbins[0:_id_])
    plt.plot(r[sum_bin:sum_bin+N],m_in_r[sum_bin:sum_bin+N],'b-',label = r'$M_{inr}$')
    plt.plot(r[sum_bin:sum_bin+N],m_star[sum_bin:sum_bin+N],'r-',label = r'$M_\ast$')
    plt.plot(r[sum_bin:sum_bin+N],m_gas[sum_bin:sum_bin+N],'g-',label = r'$M_{gas}$')
    
    m_b = m_gas[sum_bin:sum_bin+N] + m_star[sum_bin:sum_bin+N]
    plt.plot(r[sum_bin:sum_bin+N],m_b,'m-',label = r'$M_b$')
    m_h = m_in_r[sum_bin:sum_bin+N] - m_b
    plt.plot(r[sum_bin:sum_bin+N],m_h,'k-',label = r'$M_d$')
    
    plt.legend(loc=4)
    plt.xlabel(r'$r [kpc/h]$')
    plt.ylabel(r'$Mass [M_\odot /h]$')
    #plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1e0,1e3)
    plt.ylim(1e8,1e16)
    #plt.savefig('mass_profile_z_0_022.png',dpi=600)
    plt.show()
    
    #get the density
    R = r[sum_bin:sum_bin+N]
    M_in = m_in_r[sum_bin:sum_bin+N]
    M_in_gas = m_gas[sum_bin:sum_bin+N]
    rho_in = np.zeros(N-1,dtype = np.float)
    rho_h = np.zeros(N-1,dtype = np.float)
    rho_g = np.zeros(N-1,dtype = np.float)
    for k in range(N-1):
        dr = R[k+1]-R[k]
        rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
        rho_h[k] = (m_h[k+1]-m_h[k])/(4*np.pi*R[k]**2*dr)
        rho_g[k] = (M_in_gas[k+1]-M_in_gas[k])/(4*np.pi*R[k]**2*dr)
    #mean density and density
    plt.plot(R[0:-1],rho_in,'m-',label = r'$\rho_{tot}$')
    plt.plot(R[0:-1],rho_h,'k-',label = r'$\rho_d$')
    plt.plot(R[0:-1],rho_g,'g-',label = r'$\rho_{gas}$')
    plt.legend(loc=1)
    plt.xlabel(r'$r [kpc/h]$')
    plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1e0,1e3)
    plt.ylim(1e4,1e11)
    #plt.savefig('density_profile_z_0_022.png',dpi=600)
    plt.show()
    '''
    plt.plot(r[sum_bin:sum_bin+N],dens[sum_bin:sum_bin+N],'g--',label = r'$\rho$')
    plt.plot(r[sum_bin:sum_bin+N],ovdens[sum_bin:sum_bin+N],'r--',label = r'$\Delta \rho$')
    plt.legend(loc=4)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$M_\odot h^2 {kpc^3}$')
    plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    #plt.savefig('density_profile_z_0_022.png',dpi=600)
    plt.show()
    #nbins:use this to identical the particles' rows belong to the halo in the *AHF_particles file
    m_p_star = Mstar[0]/Nstar[0]
    m_p_gas = Mgas[0]/Ngas[0]
    plt.plot(r[sum_bin:sum_bin+N],n_part[sum_bin:sum_bin+N],'g',label = 'Num-particle')
    plt.legend(loc = 4)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'N')
    plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    #plt.savefig('Number_z_0_022.png',dpi=600)
    plt.show()

    plt.plot(r[sum_bin:sum_bin+N],n_part[sum_bin:sum_bin+N]/(4*np.pi*r[sum_bin:sum_bin+N]**3/3),'r',label = 'Num-density')
    plt.legend(loc = 4)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$N h^3/kpc^3$')
    plt.axvline(locpeak[_id_],ls='--',linewidth=0.5,color='blue')
    plt.yscale('log')
    plt.xscale('log')
    #plt.savefig('Number_density_z_0_022')
    plt.show() 
    '''
#try to count the particles of 



