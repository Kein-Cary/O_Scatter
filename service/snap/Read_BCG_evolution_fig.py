import matplotlib as mpl
mpl.use('Agg')
#from pygadgetreader import *
#import pygadgetreader as pygdr
import numpy as np
import h5py
#import pandas as pd
#import astropy.io.ascii as asc
#import scipy.stats as st
import matplotlib.pyplot as plt
#from handy import scatter as hsc
#import find
#import changds
def fig_BCG(resolution,g_size):
    resolution = resolution
    size_BCG = g_size
    # the value is same as the file:Read_BCG_evolution
    # read the radius of BCGs
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_R_%.3f.h5'%(size_BCG,resolution)) as f:
        r = np.array(f['a'])
    r_mu = np.array(r[0,:])
    r_gx = np.array(r[1,:])
    # read the redshift
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_redshift_MU.h5') as f:
        z0 = np.array(f['a'])
    z_mu = np.array(z0)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_redshift_GX.h5') as f:
        z1 = np.array(f['a'])    
    z_gx = np.array(z1)
    # read the data of MUSIC
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_m_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        m = np.array(f['a'])
    m_s = np.array(m[0,:])
    m_g = np.array(m[1,:])
    m_bh = np.array(m[2,:])
    m_b = m_s + m_g
    m_tot = m_b + m_bh
    
    eta_b = m_b/m_tot
    eta_s = m_s/m_tot
    eta_g = m_g/m_tot
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_rho_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        rho = np.array(f['a'])
    rho_s = np.array(rho[0,:])
    rho_g = np.array(rho[1,:])
    rho_b = rho_s + rho_g
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_mrho_MU_%.3f.h5'%(size_BCG,resolution)) as f:
        m_rho = np.array(f['a'])
    m_rho_s = np.array(m_rho[0,:])
    m_rho_g = np.array(m_rho[1,:])
    m_rho_b = m_rho_s + m_rho_g
    # read the data of GadgetX
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_m_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        m_gx = np.array(f['a'])
    m_s_gx = np.array(m_gx[0,:])
    m_g_gx = np.array(m_gx[1,:])
    m_bh_gx = np.array(m_gx[2,:])
    m_b_gx = m_s_gx + m_g_gx
    m_tot_gx = m_b_gx + m_bh_gx
    
    eta_b_gx = m_b_gx/m_tot_gx
    eta_s_gx = m_s_gx/m_tot_gx
    eta_g_gx = m_g_gx/m_tot_gx
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_rho_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        rho_gx = np.array(f['a'])
    rho_s_gx = np.array(rho_gx[0,:])
    rho_g_gx = np.array(rho_gx[1,:])
    rho_b_gx = rho_s_gx + rho_g_gx
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/BCG_%.0f_mrho_GX_%.3f.h5'%(size_BCG,resolution)) as f:
        m_rho_gx = np.array(f['a'])
    m_rho_s_gx = np.array(m_rho_gx[0,:])
    m_rho_g_gx = np.array(m_rho_gx[1,:])
    m_rho_b_gx = m_rho_s_gx + m_rho_g_gx   
    # figure the result
    plt.figure(figsize = (8,9))
    plt.plot(z_mu,m_s,'r-',label = r'$M_\ast MU$')
    plt.plot(z_mu,m_g,'g-',label = r'$M_g MU$')
    plt.plot(z_mu,m_b,'b-',label = r'$M_b MU$')
    plt.plot(z_gx,m_s_gx,'r--',label = r'$M_\ast GX$')
    plt.plot(z_gx,m_g_gx,'g--',label = r'$M_g GX$')
    plt.plot(z_gx,m_b_gx,'b--',label = r'$M_b GX$')
    handles, labels = plt.gca().get_legend_handles_labels() 
    #plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=2,
    #   ncol=3, mode = "expand",fontsize=10., borderaxespad=0.) 
    plt.legend(loc = 2)
    plt.xlabel(r'$ z $')
    plt.ylabel(r'$M [M_\odot /h]$')
    plt.yscale('log')
    plt.title(r'$BCG-m-z-Re_{%.3f}$'%resolution)
    plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_mass_growth.png',dpi = 600)
    plt.show()
    plt.close()
    
    plt.figure(figsize = (6,6))
    plt.plot(z_mu,eta_b,'b-',label = r'$\eta[M_b/M_{tot}] MU$')
    plt.plot(z_mu,eta_s,'r-',label = r'$\eta[M_\ast/M_{tot}] MU$')
    plt.plot(z_mu,eta_g,'g-',label = r'$\eta[M_g/M_{tot}] MU$')
    plt.plot(z_gx,eta_b_gx,'b--',label = r'$\eta[M_b/M_{tot}] GX$')
    plt.plot(z_gx,eta_s_gx,'r--',label = r'$\eta[M_\ast/M_{tot}] GX$')
    plt.plot(z_gx,eta_g_gx,'g--',label = r'$\eta[M_g/M_{tot}] GX$')
    handles, labels = plt.gca().get_legend_handles_labels() 
    #plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=4,
    #     ncol=3, mode = "expand",fontsize=7.5, borderaxespad=2.) 
    plt.legend(loc = 4)
    plt.xlabel(r'$z$')
    plt.ylabel(r'$\eta$')
    plt.title(r'$\eta-z-Re_{%.3f}$'%resolution)
    plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_mass_ratio.png',dpi = 600)
    plt.show()
    plt.close()
    
    N_color = np.int0(len(z_mu)) 
    #plt.figure(figsize = (8,9))
    f, (ax1, ax2, ax3) = plt.subplots(1,3,sharex = False,sharey=True)
    for k in range(N_color):
        if k % 8 == 0:
            ax1.plot(r_mu[:-1],rho_b[k,:],'--',color = mpl.cm.rainbow(k/N_color),
                     label = r'$z{%.3f}$'%z_mu[k])
            handles, labels = plt.gca().get_legend_handles_labels()
    ax1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    ax1.set_xlim(1e-1,1e2)
    ax1.set_ylim(1e4,1e11) 
    ax1.set_xticklabels(ax1.get_xticks(),rotation = 45,fontsize = 7.5)
    ax1.set_yticklabels(ax1.get_yticks(),rotation = 90,fontsize = 7.5)
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in'
    ax1.set_xlabel(r'$r [kpc/h]$')
    ax1.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_title(r'$BCG-\rho_b-r-MU-Re_{%.3f}$'%resolution, fontsize = 7.5)

    for k in range(N_color):
        if k % 8 == 0:
            ax2.plot(r_mu[:-1],rho_s[k,:],'--',color = mpl.cm.rainbow(k/N_color),
                     label = r'$z{%.3f}$'%z_mu[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    ax2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    ax2.set_xlim(1e-1,1e2)
    ax2.set_xticklabels(ax2.get_xticks(),rotation = 45,fontsize = 7.5)
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    ax2.set_xlabel(r'$r [kpc/h]$')
    ax2.set_xscale('log')
    #ax2.set_ylim(1e4,1e11)
    #ax2.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
    #ax2.set_yscale('log')
    ax2.set_title(r'$BCG-\rho_\ast-r-MU-Re_{%.3f}$'%resolution, fontsize = 7.5)
    #plt.savefig('BCG star profile evolution MU')
    #plt.show()
    #plt.close()
    #plt.figure(figsize = (8,9))
    for k in range(N_color):
        if k % 8 == 0:
            ax3.plot(r_mu[:-1],rho_g[k,:],'--',color = mpl.cm.rainbow(k/N_color),
                     label = r'$z{%.3f}$'%z_mu[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    ax3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    ax3.set_xlim(1e-1,1e2)
    ax3.set_xticklabels(ax3.get_xticks(),rotation = 45, fontsize = 7.5)
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    ax3.set_xlabel(r'$r [kpc/h]$')
    ax3.set_xscale('log')
    #ax3.set_ylim(1e4,1e11)
    #ax3.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
    #ax3.set_yscale('log')
    ax3.set_title(r'$BCG-\rho_g-r-MU-Re_{%.3f}$'%resolution, fontsize = 7.5)
    #plt.savefig('BCG gas profile evolution MU')
    #plt.show()
    #plt.close()
    plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
    plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_density_profile_evolution_MU.png',dpi = 600)
    plt.show()
    plt.close()

    N_color_gx = np.int0(len(z_gx))
    #plt.figure(figsize = (8,9))
    f, (bx1, bx2, bx3) = plt.subplots(1,3,sharex = False,sharey=True)
    for k in range(N_color_gx):
        if k % 8 == 0:
            bx1.plot(r_gx[:-1],rho_b_gx[k,:],'--',color = mpl.cm.rainbow(k/N_color_gx),
                     label = r'$z{%.3f}$'%z_gx[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    bx1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    bx1.set_xlabel(r'$r [kpc/h]$')
    bx1.set_xlim(1e-1,1e2)
    bx1.set_ylim(1e4,1e11)
    bx1.set_xticklabels(bx1.get_xticks(),rotation = 45,fontsize = 7.5)
    bx1.set_yticklabels(bx1.get_yticks(),rotation = 90,fontsize = 7.5) 
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    bx1.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
    bx1.set_yscale('log')
    bx1.set_xscale('log')
    bx1.set_title(r'$BCG-\rho_b-r-GX-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG baryon profile evolution GX')
    #plt.show()
    #plt.close()
    #plt.figure(figsize = (8,9))
    for k in range(N_color_gx):
        if k % 8 == 0:
            bx2.plot(r_gx[:-1],rho_s_gx[k,:],'--',color = mpl.cm.rainbow(k/N_color_gx),
                     label = r'$z{%.3f}$'%z_gx[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    bx2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    bx2.set_xlim(1e-1,1e2)
    bx2.set_xticklabels(bx2.get_xticks(),rotation = 45,fontsize = 7.5)   
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    bx2.set_xlabel(r'$r [kpc/h]$')
    bx2.set_xscale('log')
    #bx2.set_ylim(1e4,1e11)    
    #bx2.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
    #bx2.set_yscale('log')
    bx2.set_title(r'$BCG-\rho_\ast-r-GX-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG star profile evolution GX')
    #plt.show()
    #plt.close()
    #plt.figure(figsize = (8,9))
    for k in range(N_color_gx):
        if k % 8 == 0:
            bx3.plot(r_gx[:-1],rho_g_gx[k,:],'--',color = mpl.cm.rainbow(k/N_color_gx),
                     label = r'$z{%.3f}$'%z_gx[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    bx3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    bx3.set_xlim(1e-1,1e2)
    bx3.set_xticklabels(bx3.get_xticks(),rotation = 45,fontsize = 7.5)   
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    bx3.set_xlabel(r'$r [kpc/h]$')
    bx3.set_xscale('log')
    #bx3.set_ylim(1e4,1e11)
    #bx3.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
    #bx3.set_yscale('log')
    bx3.set_title(r'$BCG-\rho_g-r-GX-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG gas profile evolution GX')
    #plt.show()
    #plt.close()
    plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
    plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_density_profile_evolution_GX.png',dpi = 600)
    plt.show()
    plt.close()
#################
    N_color = np.int0(len(z_mu)) 
    #plt.figure(figsize = (8,9))
    f, (cx1, cx2, cx3) = plt.subplots(1,3,sharex = False,sharey=True)
    for k in range(N_color):
        if k % 8 == 0:
            cx1.plot(r_mu,m_rho_b[k,:],'--',color = mpl.cm.rainbow(k/N_color),
                     label = r'$z{%.3f}$'%z_mu[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    cx1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    cx1.set_xlim(1e-1,1e2)
    cx1.set_ylim(1e4,1e11)
    cx1.set_xticklabels(cx1.get_xticks(),rotation = 45,fontsize = 7.5)
    cx1.set_yticklabels(cx1.get_yticks(),rotation = 90,fontsize = 7.5) 
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    cx1.set_xlabel(r'$r [kpc/h]$')
    cx1.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
    cx1.set_yscale('log')
    cx1.set_xscale('log')
    cx1.set_title(r'$BCG-\bar{\rho_b}-r-MU-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG mean baryon profile evolution MU')
    #plt.show()
    #plt.close()
    #plt.figure(figsize = (8,9))
    for k in range(N_color):
        if k % 8 == 0:
            cx2.plot(r_mu,m_rho_s[k,:],'--',color = mpl.cm.rainbow(k/N_color),
                     label = r'$z{%.3f}$'%z_mu[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    cx2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    cx2.set_xlim(1e-1,1e2)
    cx2.set_xticklabels(cx2.get_xticks(),rotation = 45,fontsize = 7.5)
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    cx2.set_xlabel(r'$r [kpc/h]$')
    cx2.set_xscale('log')
    #cx2.set_ylim(1e4,1e11)
    #cx2.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
    #cx2.set_yscale('log')
    cx2.set_title(r'$BCG-\bar{\rho_\ast}-r-MU-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG mean star profile evolution MU')
    #plt.show()
    #plt.close()
    #plt.figure(figsize = (8,9))
    for k in range(N_color):
        if k % 8 == 0:
            cx3.plot(r_mu,m_rho_g[k,:],'--',color = mpl.cm.rainbow(k/N_color),
                     label = r'$z{%.3f}$'%z_mu[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    cx3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    cx3.set_xlim(1e-1,1e2)
    cx3.set_xticklabels(cx3.get_xticks(),rotation = 45,fontsize = 7.5)
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    cx3.set_xlabel(r'$r [kpc/h]$')
    cx3.set_xscale('log')
    #cx3.set_ylim(1e4,1e11)
    #cx3.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
    #cx3.set_yscale('log')
    cx3.set_title(r'$BCG-\bar{\rho_g}-r-MU-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG mean gas profile evolution MU')
    #plt.show()
    #plt.close()
    plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
    plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_mean_density_profile_evolution_MU.png',dpi = 600)
    plt.show()
    plt.close()
    
    N_color_gx = np.int0(len(z_gx))
    #plt.figure(figsize = (8,9))
    f, (dx1, dx2, dx3) = plt.subplots(1,3,sharex = False,sharey=True)
    for k in range(N_color_gx):
        if k % 8 == 0:
            dx1.plot(r_gx,m_rho_b_gx[k,:],'--',color = mpl.cm.rainbow(k/N_color_gx),
                     label = r'$z{%.3f}$'%z_gx[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    dx1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad=0.)
    dx1.set_xlabel(r'$r [kpc/h]$')
    dx1.set_xlim(1e-1,1e2)
    dx1.set_ylim(1e4,1e11)
    dx1.set_xticklabels(dx1.get_xticks(),rotation = 45,fontsize = 7.5)
    dx1.set_yticklabels(dx1.get_yticks(),rotation = 90,fontsize = 7.5) 
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    dx1.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
    dx1.set_yscale('log')
    dx1.set_xscale('log')
    dx1.set_title(r'$BCG-\bar{\rho_b}-r-GX-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG mean baryon profile evolution GX')
    #plt.show()
    #plt.close()
    #plt.figure(figsize = (8,9))
    for k in range(N_color_gx):
        if k % 8 == 0:
            dx2.plot(r_gx,m_rho_s_gx[k,:],'--',color = mpl.cm.rainbow(k/N_color_gx),
                     label = r'$z{%.3f}$'%z_gx[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    dx2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize= 5., borderaxespad= 0.)
    dx2.set_xlim(1e-1,1e2)
    dx2.set_xticklabels(dx2.get_xticks(),rotation = 45,fontsize = 7.5)
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in'     
    dx2.set_xlabel(r'$r [kpc/h]$')
    dx2.set_xscale('log')
    #dx2.set_ylim(1e4,1e11)    
    #dx2.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
    #dx2.set_yscale('log')
    dx2.set_title(r'$BCG-\bar{\rho_\ast}-r-GX-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG mean star profile evolution GX')
    #plt.show()
    #plt.close()
    #plt.figure(figsize = (8,9))
    for k in range(N_color_gx):
        if k % 8 == 0:
            dx3.plot(r_gx,m_rho_g_gx[k,:],'--',color = mpl.cm.rainbow(k/N_color_gx),
                     label = r'$\bar{\rho_g}-{%.3f}$'%z_gx[k])
            handles, labels = plt.gca().get_legend_handles_labels() 
    dx3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
               ncol= 3, mode="expand", fontsize = 5., borderaxespad=0.)
    dx3.set_xlim(1e-1,1e2)
    dx3.set_xticklabels(dx3.get_xticks(),rotation = 45,fontsize = 7.5)   
    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    dx3.set_xlabel(r'$r [kpc/h]$')
    dx3.set_xscale('log')
    #dx3.set_ylim(1e4,1e11)
    #dx3.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
    #dx3.set_yscale('log')
    dx3.set_title(r'$BCG-\bar{\rho_g}-r-GX-Re_{%.3f}$'%resolution,fontsize = 7.5)
    #plt.savefig('BCG mean gas profile evolution GX')
    #plt.show()
    #plt.close()
    plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
    plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_mean_density_profile_evolution_GX.png',dpi = 600)
    plt.show()
    plt.close()
    return
#fig_BCG(resolution = True,trace = True)
