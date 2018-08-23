#try to figure the profile of GX simulation
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
####part 0:read the redshift,to show the result in different z
with h5py.File('Redshift_GX.h5') as f:
    y0 = f['a']
    z = np.array(y0)
with h5py.File('main_tree_GX.h5') as f:
    y1 = f['a']
    main_tree = np.array(y1)
##assuming now get the first one halo at z=0,and see the profile at z
goal_halo = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
goalhalo = np.array(goal_halo['col1'])
goal_id = np.int64(goalhalo[0])##set the goal halo
#set the range of z,during this range to see the profile evolution
pp = z<=1.0
z_rang = z[pp]
l_rang = len(z_rang)
plt.figure(figsize=(8,9))
#plt.figure(figsize=(8,8))#for density
for q in range(l_rang):
    if q % 6 == 0:
        goal_z = z_rang[q]
        IX = z == goal_z
        IY = IX.tolist()
        No_site_z = IY.index(True)
        ia = main_tree[:,0] == goal_id
        ib = ia.tolist()
        No_site_id = ib.index(True)
        #get the goal main progenitor
        goal_progenitor = main_tree[No_site_id,No_site_z]
        ###############part1
        M_halo = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos'%z[No_site_z],
                          converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
                          ## this line try to set the datatype of the first coloum
        # try to get the base property of halo 128000000000001
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
        ID = np.int64(ID) 
        #site the goal ID at z!=0 from main_tree
        ip = ID == goal_progenitor
        iq = ip.tolist()
        _ix = iq.index(True)
        ix = Host == ID[_ix]################part2
        h_profile = pd.read_table('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_profiles'%z[No_site_z],dtype = np.float)
        r = np.array(h_profile['#r(1)'])
        r = np.abs(r)#to make sure that r is positive
        n_part = np.array(h_profile['npart(2)'])
        m_in_r = np.array(h_profile['M_in_r(3)'])
        m_star = np.array(h_profile['M_star(26)'])
        m_gas = np.array(h_profile['M_gas(25)'])
        dens = np.array(h_profile['dens(5)'])
        ovdens = np.array(h_profile['ovdens(4)'])
        dm = m_in_r - m_star - m_gas
        #figure
        N = np.int(Nbins[_ix])
        if _ix==0 :
            m_b = m_gas[0:N] + m_star[0:N]
            m_h = m_in_r[0:N] - m_b
            '''
            #plt.plot(r[0:N],m_in_r[0:N],':',color = mpl.cm.rainbow(q/l_rang),label = r'$M_{inr} -%.3f$'%goal_z)
            #plt.plot(r[0:N],m_star[0:N],'--',color = mpl.cm.rainbow(q/l_rang),label = r'$M_\ast -%.3f$'%goal_z)
            plt.plot(r[0:N],m_gas[0:N],'.-',color = mpl.cm.rainbow(q/l_rang),label = r'$M_{gas} -%.3f$'%goal_z)
            #plt.plot(r[0:N],m_b,'-',color = mpl.cm.rainbow(q/l_rang),label = r'$M_b -%.3f$'%goal_z)
            plt.plot(r[0:N],m_h,':.',color = mpl.cm.rainbow(q/l_rang),label = r'$M_d -%.3f$'%goal_z)
            '''
            #get the density
            R = r[0:N]
            M_in = m_in_r[0:N]
            M_in_gas = m_gas[0:N]
            M_star = m_star[0:N]
            M_b = M_in_gas+M_star
            rho_in = np.zeros(N-1,dtype = np.float)
            rho_h = np.zeros(N-1,dtype = np.float)
            rho_g = np.zeros(N-1,dtype = np.float)
            rho_s = np.zeros(N-1,dtype = np.float)
            rho_b = np.zeros(N-1,dtype = np.float)
            for k in range(N-1):
                dr = R[k+1]-R[k]
                rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
                rho_h[k] = (m_h[k+1]-m_h[k])/(4*np.pi*R[k]**2*dr)
                rho_g[k] = (M_in_gas[k+1]-M_in_gas[k])/(4*np.pi*R[k]**2*dr)
                rho_s[k] = (M_star[k+1]-M_star[k])/(4*np.pi*R[k]**2*dr)
                rho_b[k] = (M_b[k+1]-M_b[k])/(4*np.pi*R[k]**2*dr)
            #plt.plot(R[0:-1],rho_in, ':',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_{tot} -%.3f$'%goal_z)
            plt.plot(R[0:-1],rho_h, ':.',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_d -%.3f$'%goal_z)
            #plt.plot(R[0:-1],rho_g, '.-',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_g -%.3f$'%goal_z)
            plt.plot(R[0:-1],rho_s, '--',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_\ast -%.3f$'%goal_z)
            #plt.plot(R[0:-1],rho_b, '-',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_b -%.3f$'%goal_z)
        else:
            sum_bin = np.sum(Nbins[0:_ix])
            m_b = m_gas[sum_bin:sum_bin+N] + m_star[sum_bin:sum_bin+N]
            m_h = m_in_r[sum_bin:sum_bin+N] - m_b
            '''
            #plt.plot(r[sum_bin:sum_bin+N],m_in_r[sum_bin:sum_bin+N],':',color = mpl.cm.rainbow(q/l_rang),label = r'$M_{inr} -%.3f$'%goal_z)
            #plt.plot(r[sum_bin:sum_bin+N],m_star[sum_bin:sum_bin+N],'--',color = mpl.cm.rainbow(q/l_rang),label = r'$M_\ast -%.3f$'%goal_z)
            plt.plot(r[sum_bin:sum_bin+N],m_gas[sum_bin:sum_bin+N],'.-',color = mpl.cm.rainbow(q/l_rang),label = r'$M_{gas} -%.3f$'%goal_z)
            #plt.plot(r[sum_bin:sum_bin+N],m_b,'-',color = mpl.cm.rainbow(q/l_rang),label = r'$M_b -%.3f$'%goal_z)
            plt.plot(r[sum_bin:sum_bin+N],m_h,':.',color = mpl.cm.rainbow(q/l_rang),label = r'$M_d -%.3f$'%goal_z)
            '''
            #get the density
            R = r[sum_bin:sum_bin+N]
            M_in = m_in_r[sum_bin:sum_bin+N]
            M_in_gas = m_gas[sum_bin:sum_bin+N]
            M_star = m_star[sum_bin:sum_bin+N]
            M_b = M_in_gas+M_star
            rho_in = np.zeros(N-1,dtype = np.float)
            rho_h = np.zeros(N-1,dtype = np.float)
            rho_g = np.zeros(N-1,dtype = np.float)
            rho_s = np.zeros(N-1,dtype = np.float)
            rho_b = np.zeros(N-1,dtype = np.float)
            for k in range(N-1):
                dr = R[k+1]-R[k]
                rho_in[k] = (M_in[k+1]-M_in[k])/(4*np.pi*R[k]**2*dr)
                rho_h[k] = (m_h[k+1]-m_h[k])/(4*np.pi*R[k]**2*dr)
                rho_g[k] = (M_in_gas[k+1]-M_in_gas[k])/(4*np.pi*R[k]**2*dr)
                rho_s[k] = (M_star[k+1]-M_star[k])/(4*np.pi*R[k]**2*dr)
                rho_b[k] = (M_b[k+1]-M_b[k])/(4*np.pi*R[k]**2*dr)
            #plt.plot(R[0:-1],rho_in,':',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_{tot} -%.3f$'%goal_z)
            plt.plot(R[0:-1],rho_h,':.',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_d -%.3f$'%goal_z)
            #plt.plot(R[0:-1],rho_g,'.-',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_g -%.3f$'%goal_z)
            plt.plot(R[0:-1],rho_s,'--',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_\ast -%.3f$'%goal_z)
            #plt.plot(R[0:-1],rho_b,'-',color = mpl.cm.rainbow(q/l_rang),label = r'$\rho_b -%.3f$'%goal_z)
'''            
plt.ylabel(r'$M [M_\odot /h]$')        
plt.xlabel(r'$r [kpc/h]$')
#plt.legend(loc=4)
plt.legend(bbox_to_anchor=(-0.05, 1., 1.15, 0.), loc=4,
           ncol=6, mode="expand", fontsize=7.5, borderaxespad=2.)
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e0,1e3)
plt.ylim(1e8,1e16)
plt.savefig('mass evolution GX',dpi = 600)
plt.show()
'''
plt.ylabel(r'$\rho [M_\odot h^2 {kpc^3}]$')
#plt.legend(loc=1)
plt.legend(bbox_to_anchor=(-0.05, 1., 1.15, 0.), loc=4,
           ncol=6, mode="expand", fontsize=7.5, borderaxespad=2.)
plt.xlabel(r'$r [kpc/h]$')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e0,1e3)
plt.ylim(1e4,1e11)
plt.savefig('density evolution GX',dpi = 600)
plt.show()
