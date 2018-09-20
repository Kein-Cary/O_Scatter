# this file try to figure the mass evolution of halo and BCG
import matplotlib as mpl
mpl.use('Agg')
from pygadgetreader import *
import pygadgetreader as pygdr
import numpy as np
import h5py
import pandas as pd
import astropy.io.ascii as asc
import scipy.stats as st
import matplotlib.pyplot as plt
from handy import scatter as hsc
import find
def compare_fig(resolution,g_size):
    resolution = resolution
    size_BCG = g_size
    # read the halo data
    with h5py.File('/home/cxkttwl/Scatter_data_read/total_mass_MUSIC.h5') as f:
        y0 = f['a']
        mass_MUSIC = np.array(y0)
    with h5py.File('/home/cxkttwl/Scatter_data_read/total_mass_GX.h5') as f:
        y1 = f['a']
        mass_GX = np.array(y1)
    with h5py.File('/home/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/redshift.h5') as f:
        y2 = f['a']
        redshift = np.array(y2)
    with h5py.File('/home/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/Redshift_GX.h5') as f:
        y3 = f['a']
        redshift_gx = np.array(y3)
    #data of MUSIC simulation
    halo_mass = np.array(mass_MUSIC[0,:])
    star_mass = np.array(mass_MUSIC[1,:])
    gas_mass = np.array(mass_MUSIC[2,:])
    baryon_m = np.array(mass_MUSIC[3,:])
    dark_mass = np.array(mass_MUSIC[4,:])
    #data of GadgetX simulation
    halo_mass_gx = np.array(mass_GX[0,:])
    star_mass_gx = np.array(mass_GX[1,:])
    gas_mass_gx = np.array(mass_GX[2,:])
    baryon_m_gx = np.array(mass_GX[3,:])
    dark_mass_gx = np.array(mass_GX[4,:])
    #figure
    ia = redshift<=1.0 
    z_MUSIC = np.array(redshift[ia])
    ia_gx = redshift_gx<=1.0 
    z_GX = np.array(redshift_gx[ia_gx])
    '''
    #find the mass change in a short time
    N = len(z_MUSIC[0:-1])
    M = len(z_GX[0:-1])
    deta_mh = np.zeros(N,dtype = np.float)
    deta_ms = np.zeros(N,dtype = np.float)
    deta_mg = np.zeros(N,dtype = np.float)
    deta_mb = np.zeros(N,dtype = np.float)
    deta_md = np.zeros(N,dtype = np.float)
    #for GadgetX 
    deta_mh_gx = np.zeros(M,dtype = np.float)
    deta_ms_gx = np.zeros(M,dtype = np.float)
    deta_mg_gx = np.zeros(M,dtype = np.float)
    deta_mb_gx = np.zeros(M,dtype = np.float)
    deta_md_gx = np.zeros(M,dtype = np.float)
    for k in range(N):
        dz = z_MUSIC[k+1]-z_MUSIC[k]
        deta_mh[k] = (halo_mass[k]-halo_mass[k+1])/dz
        deta_ms[k] = (star_mass[k]-star_mass[k+1])/dz
        deta_mg[k] = (gas_mass[k]-gas_mass[k+1])/dz
        deta_mb[k] = (baryon_m[k]-baryon_m[k+1])/dz
        deta_md[k] = (dark_mass[k]-dark_mass[k+1])/dz
    for k in range(M):
        dz_gx = z_GX[k+1]-z_GX[k]
        deta_mh_gx[k] = (halo_mass_gx[k]-halo_mass_gx[k+1])/dz_gx
        deta_ms_gx[k] = (star_mass_gx[k]-star_mass_gx[k+1])/dz_gx
        deta_mg_gx[k] = (gas_mass_gx[k]-gas_mass_gx[k+1])/dz_gx
        deta_mb_gx[k] = (baryon_m_gx[k]-baryon_m_gx[k+1])/dz_gx
        deta_md_gx[k] = (dark_mass_gx[k]-dark_mass_gx[k+1])/dz_gx
    '''
    ratio_gas_b = gas_mass/baryon_m
    ratio_gas_d = gas_mass/dark_mass
    ratio_star_d = star_mass/dark_mass
    ratio_star_b = star_mass/baryon_m
    ratio_b = baryon_m/halo_mass
    ratio_d = dark_mass/halo_mass
    ratio_gas_b_gx = gas_mass_gx/baryon_m_gx
    ratio_gas_d_gx = gas_mass_gx/dark_mass_gx
    ratio_star_d_gx = star_mass_gx/dark_mass_gx
    ratio_star_b_gx = star_mass_gx/baryon_m_gx
    ratio_b_gx = baryon_m_gx/halo_mass_gx
    ratio_d_gx = dark_mass_gx/halo_mass_gx
    ######################################
    # read the BCG data
    with h5py.File('BCG_%.0f_R_%.3f.h5'%(size_BCG,resolution)) as f:
        r = np.array(f['a'])
    r_mu = np.array(r[0,:])
    r_gx = np.array(r[1,:])
    # read the redshift
    with h5py.File('BCG_redshift.h5') as f:
        z = np.array(f['a'])
    z_mu = np.array(z[0,:])
    z_gx = np.array(z[1,:])
    # read the data of MUSIC
    with h5py.File('BCG_%.0f_m_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        m = np.array(f['a'])
    m_s = np.array(m[0,:])
    m_g = np.array(m[1,:])
    m_bh = np.array(m[2,:])
    m_b = m_s + m_g
    m_tot = m_b + m_bh
    eta_b = m_b/m_tot
    eta_s = m_s/m_tot
    eta_g = m_g/m_tot
    # read the data of GadgetX
    with h5py.File('BCG_%.0f_m_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        m_gx = np.array(f['a'])
    m_s_gx = np.array(m_gx[0,:])
    m_g_gx = np.array(m_gx[1,:])
    m_bh_gx = np.array(m_gx[2,:])
    m_b_gx = m_s_gx + m_g_gx
    m_tot_gx = m_b_gx + m_bh_gx  
    eta_b_gx = m_b_gx/m_tot_gx
    eta_s_gx = m_s_gx/m_tot_gx
    eta_g_gx = m_g_gx/m_tot_gx
    '''
    N_bcg = len(z_mu[0:-1])
    M_bcg = len(z_gx[0:-1])
    deta_s = np.zeros(N_bcg,dtype = np.float)
    deta_g = np.zeros(N_bcg,dtype = np.float)
    deta_b = np.zeros(N_bcg,dtype = np.float)
    #for GadgetX 
    deta_s_gx = np.zeros(M_bcg,dtype = np.float)
    deta_g_gx = np.zeros(M_bcg,dtype = np.float)
    deta_b_gx = np.zeros(M_bcg,dtype = np.float)
    for k in range(N_bcg):
        dz_bcg = z_mu[k+1]-z_mu[k]
        deta_s[k] = (m_s[k]-m_s[k+1])/dz_bcg
        deta_g[k] = (m_g[k]-m_g[k+1])/dz_bcg
        deta_b[k] = (m_b[k]-m_b[k+1])/dz_bcg
    for k in range(M_bcg):
        dz_bcg_gx = z_gx[k+1]-z_gx[k]
        deta_s_gx[k] = (m_s_gx[k]-m_s_gx[k+1])/dz_bcg_gx
        deta_g_gx[k] = (m_g_gx[k]-m_g_gx[k+1])/dz_bcg_gx
        deta_b_gx[k] = (m_b_gx[k]-m_b_gx[k+1])/dz_bcg_gx
    ## fig_out the result
    plt.figure(figsize=(8,9))
    plt.plot(z_MUSIC[0:-1],deta_ms,'r-',label = r'$dM_\ast/dz[MU]$')
    plt.plot(z_MUSIC[0:-1],deta_mg,'g-',label = r'$dM_g/dz[MU]$')
    plt.plot(z_MUSIC[0:-1],deta_mb,'m-',label = r'$dM_b/dz[MU]$')
    plt.plot(z_MUSIC[0:-1],deta_md,'k-',label = r'$dM_d/dz[MU]$')   
    plt.plot(z_GX[0:-1],deta_ms_gx,'r--',label = r'$dM_\ast/dz[GX]$')
    plt.plot(z_GX[0:-1],deta_mg_gx,'g--',label = r'$dM_g/dz[GX]$')
    plt.plot(z_GX[0:-1],deta_mb_gx,'m--',label = r'$dM_b/dz[GX]$')
    plt.plot(z_GX[0:-1],deta_md_gx,'k--',label = r'$dM_d/dz[GX]$')
    plt.plot(z_mu[0:-1],deta_s,'r.:',label = r'$M_\ast[MU-bcg]$')
    plt.plot(z_mu[0:-1],deta_g,'g.:',label = r'$M_g[MU-bcg]$')
    plt.plot(z_mu[0:-1],deta_b,'b.:',label = r'$M_b[MU-bcg]$')
    plt.plot(z_gx[0:-1],deta_s_gx,'r:',label = r'$M_\ast[GX-bcg]$')
    plt.plot(z_gx[0:-1],deta_g_gx,'g:',label = r'$M_g[GX-bcg]$')
    plt.plot(z_gx[0:-1],deta_b_gx,'b:',label = r'$M_b[GX-bcg]$')
    handles, labels = plt.gca().get_legend_handles_labels() 
    plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=4,
       ncol= 4, mode = "expand",fontsize=10., borderaxespad=0.) 
    plt.yscale('log')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$M[M_\odot /h]/z$')
    plt.title(r'$BCG-M_{raise}-z$','right')
    plt.savefig('M growth',dpi = 600)
    plt.show()
    plt.close()
    '''
    plt.figure(figsize=(8,9))
    plt.plot(z_MUSIC,halo_mass,'b-',label = r'$M_h[MU]$')
    plt.plot(z_MUSIC,star_mass,'r-',label = r'$M_s[MU]$')
    plt.plot(z_MUSIC,gas_mass,'g-',label = r'$M_g[MU]$')
    plt.plot(z_MUSIC,baryon_m,'m-',label = r'$M_b[MU]$')
    plt.plot(z_MUSIC,dark_mass,'k-',label = r'$M_d[MU]$')
    plt.plot(z_GX,halo_mass_gx,'b--',label = r'$M_h[GX]$')
    plt.plot(z_GX,star_mass_gx,'r--',label = r'$M_s[GX]$')
    plt.plot(z_GX,gas_mass_gx,'g--',label = r'$M_g[GX]$')
    plt.plot(z_GX,baryon_m_gx,'m--',label = r'$M_b[GX]$')
    plt.plot(z_GX,dark_mass_gx,'k--',label = r'$M_d[GX]$')
    plt.plot(z_mu,m_s,'r.-',label = r'$M_\ast[MU-bcg]$')
    plt.plot(z_mu,m_g,'g.-',label = r'$M_g[MU-bcg]$')
    plt.plot(z_mu,m_b,'b.-',label = r'$M_b[MU-bcg]$')
    plt.plot(z_gx,m_s_gx,'r:',label = r'$M_\ast[GX-bcg]$')
    plt.plot(z_gx,m_g_gx,'g:',label = r'$M_g[GX-bcg]$')
    plt.plot(z_gx,m_b_gx,'b:',label = r'$M_b[GX-bcg]$')
    handles, labels = plt.gca().get_legend_handles_labels() 
    plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=4,
       ncol= 4, mode = "expand",fontsize=10., borderaxespad=0.) 
    plt.yscale('log')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$M[M_\odot /h]$')
    #plt.title(r'$BCG-M-z$',loc = 'right')
    plt.text(0.5,3e15,'BCG&halo-M-Z',fontsize=20.)
    plt.savefig('M evolution',dpi = 600)
    plt.show()
    plt.close()
    '''
    plt.figure(figsize=(10,10))
    plt.plot(z_MUSIC,ratio_gas_b,'b-',label = r'$\eta[M_g-M_b]-MU$')
    plt.plot(z_MUSIC,ratio_d,'k-',label = r'$\eta[M_d-M_h]-MU$')
    plt.plot(z_GX,ratio_gas_b_gx,'b--',label = r'$\eta[M_g-M_b]-GX$')
    plt.plot(z_GX,ratio_d_gx,'k--',label = r'$\eta[M_d-M_h]-GX$')
    plt.plot(z_MUSIC,ratio_gas_d,'g-',label = r'$\eta[M_g-M_d]-MU$')
    plt.plot(z_MUSIC,ratio_star_d,'r-',label = r'$\eta[M_\ast-M_d]-MU$')
    plt.plot(z_MUSIC,ratio_star_b,'m-',label = r'$\eta[M_\ast-M_b]-MU$')
    plt.plot(z_MUSIC,ratio_b,'y-',label = r'$\eta[M_b-M_h]-MU$')
    plt.plot(z_GX,ratio_gas_d_gx,'g--',label = r'$\eta[M_g-M_d]-GX$')
    plt.plot(z_GX,ratio_star_d_gx,'r--',label = r'$\eta[M_\ast-M_d]-GX$')
    plt.plot(z_GX,ratio_star_b_gx,'m--',label = r'$\eta[M_\ast-M_b]-GX$')
    plt.plot(z_GX,ratio_b_gx,'y--',label = r'$\eta[M_b-M_h]-GX$')
    plt.plot(z_mu,eta_s,'r.-',label = r'$\eta[M_\ast-M_{tot}]-MU_{bcg}$')
    plt.plot(z_mu,eta_g,'g.-',label = r'$\eta[M_g-M_{tot}]-MU_{bcg}$')
    plt.plot(z_mu,eta_b,'b.-',label = r'$\eta[M_b-M_{tot}]-MU_{bcg}$')
    plt.plot(z_gx,eta_s_gx,'r:',label = r'$\eta[M_\ast-M_{tot}]-GX_{bcg}$')
    plt.plot(z_gx,eta_g_gx,'g:',label = r'$\eta[M_g-M_{tot}]-GX_{bcg}$')
    plt.plot(z_gx,eta_b_gx,'b:',label = r'$\eta[M_b-M_{tot}]-GX_{bcg}$') 
    handles, labels = plt.gca().get_legend_handles_labels() 
    plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=4,
       ncol=5, mode = "expand",fontsize=10., borderaxespad=0.)
    plt.xlabel(r'$z$')
    plt.ylabel(r'$\eta_m$')
    plt.title(r'$BCG-\eta_m-z$',loc = 'right')
    plt.savefig('Ratio variables',dpi = 600)
    plt.show()
    plt.close()
    '''
    f,(ax1,ax2,ax3) = plt.subplots(3, sharex=True, sharey=False)
    ax1.plot(z_MUSIC,ratio_gas_b,'b-',label = r'$\eta -M_g-M_b$')
    ax1.plot(z_MUSIC,ratio_d,'k-',label = r'$\eta -M_d-M_h$')
    ax1.plot(z_GX,ratio_gas_b_gx,'b--',label = r'$\eta -M_g-M_b -GX$')
    ax1.plot(z_GX,ratio_d_gx,'k--',label = r'$\eta -M_d-M_h -GX$')
    handles, labels = plt.gca().get_legend_handles_labels()
    ax1.legend(bbox_to_anchor=(0., 0.85, 1., 0.), loc=2,
         ncol=2,fontsize=7.5, borderaxespad=0.) 
    plt.ylabel(r'$\eta$')
    plt.ylim(0.75,0.95)
    ax1.set_title(r'$\eta - z $')
    ax2.plot(z_MUSIC,ratio_gas_d,'g-',label = r'$\eta[M_g-M_d]$')
    ax2.plot(z_MUSIC,ratio_star_d,'r-',label = r'$\eta[M_\ast-M_d]$')
    ax2.plot(z_MUSIC,ratio_star_b,'m-',label = r'$\eta[M_\ast-M_b]$')
    ax2.plot(z_MUSIC,ratio_b,'y-',label = r'$\eta[M_b-M_h]$')
    ax2.plot(z_GX,ratio_gas_d_gx,'g--',label = r'$\eta[M_g-M_d]GX$')
    ax2.plot(z_GX,ratio_star_d_gx,'r--',label = r'$\eta[M_\ast-M_d]GX$')
    ax2.plot(z_GX,ratio_star_b_gx,'m--',label = r'$\eta[M_\ast-M_b]GX$')
    ax2.plot(z_GX,ratio_b_gx,'y--',label = r'$\eta[M_b-M_h]GX$')
    handles, labels = plt.gca().get_legend_handles_labels()
    ax2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=3,
         ncol=4,mode = "expand",fontsize=7.5, borderaxespad=0.) 
    plt.ylim(0.,0.5)
    plt.ylabel(r'$\eta$')
    plt.xlabel(r'$ z $')
    f.subplots_adjust(hspace=0.35)
    plt.setp([a.get_xticklabels() for a in f.axes[:-2]], visible=False)
    
    ax3.plot(z_mu,eta_s,'r.-',label = r'$\eta[M_\ast-M_{tot}]-MU_{bcg}$')
    ax3.plot(z_mu,eta_g,'g.-',label = r'$\eta[M_g-M_{tot}]-MU_{bcg}$')
    ax3.plot(z_mu,eta_b,'b.-',label = r'$\eta[M_b-M_{tot}]-MU_{bcg}$')
    ax3.plot(z_gx,eta_s_gx,'r:',label = r'$\eta[M_\ast-M_{tot}]-GX_{bcg}$')
    ax3.plot(z_gx,eta_g_gx,'g:',label = r'$\eta[M_g-M_{tot}]-GX_{bcg}$')
    ax3.plot(z_gx,eta_b_gx,'b:',label = r'$\eta[M_b-M_{tot}]-GX_{bcg}$') 
    handles, labels = plt.gca().get_legend_handles_labels() 
    ax3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=3,
         ncol=3,mode = "expand",fontsize=7.5, borderaxespad=0.) 
    plt.ylim(0.,1.1)
    plt.ylabel(r'$\eta$')
    plt.xlabel(r'$ z $')
    plt.savefig('ratio Mass',dpi = 600)
    plt.show()
    plt.close()
    return r_mu,r_gx
#compare_fig(resolution = True,g_size = True)