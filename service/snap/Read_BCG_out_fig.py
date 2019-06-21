# this file try to figure the BCG in a given halo
import matplotlib as mpl
mpl.use('Agg')
#from pygadgetreader import *
import pygadgetreader as pygdr
import numpy as np
import h5py
import astropy.io.ascii as asc
import matplotlib.pyplot as plt
#from handy import scatter as hsc
import find
import changds 
def find_BCG_fig(ip,resolution,g_size):
    resolution = resolution
    size_BCG = g_size
    _id_ = ip
    Nr = 1./resolution
    # to make sure the bins density is the same
    N_size = 50.*(size_BCG/100.0)
    Nsize = np.int0(np.ceil(N_size))
    R_BCG = size_BCG
    r_bcg = np.logspace(1e-3,np.log10(R_BCG),Nsize)
    R_BCG_gx = size_BCG
    r_bcg_gx = np.logspace(1e-3,np.log10(R_BCG_gx),Nsize)
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_z_mu.h5') as f:
        zmu = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_z_gx.h5') as f:
        zgx = np.array(f['a'])  
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_s_mu.h5') as f:
        m_s = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_g_mu.h5') as f:
        m_g = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_b_mu.h5') as f:
        m_b = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_s_mu.h5') as f:
        rho_s = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_g_mu.h5') as f:
        rho_g = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_b_mu.h5') as f:
        rho_b = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_s_mu.h5') as f:
        m_rho_s = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_g_mu.h5') as f:
        m_rho_g = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_b_mu.h5') as f:
        m_rho_b = np.array(f['a'])
    ###############################
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_s_gx.h5') as f:
        m_s_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_g_gx.h5') as f:
        m_g_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_m_b_gx.h5') as f:
        m_b_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_s_gx.h5') as f:
        rho_s_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_g_gx.h5') as f:
        rho_g_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_rho_b_gx.h5') as f:
        rho_b_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_s_gx.h5') as f:
        m_rho_s_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_g_gx.h5') as f:
        m_rho_g_gx = np.array(f['a'])
    with h5py.File(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/h5_data/Read_out_mrho_b_gx.h5') as f:
        m_rho_b_gx = np.array(f['a'])
    #read the meger tree to get the star time of each halo
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
    #figure the MUSIC simulation
    for k in range(len(zmu)):
        plt.figure()
        f, (ax1, ax2, ax3) = plt.subplots(1,3,sharex = False,sharey=True)
        ax1.plot(r_bcg, m_g[k,:], 'g-', label=r'$M_g$')
        ax1.legend(loc = 2)
        ax1.set_ylim(1e8,1e13)
        ax1.set_xlim(1e-1, np.max(r_bcg))
        ax1.set_xticklabels(ax1.get_xticks(),rotation = 45,fontsize = 7.5)
        ax1.set_yticklabels(ax1.get_yticks(),rotation = 90,fontsize = 7.5)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.set_xlabel(r'$r-kpc/h$')
        ax1.set_ylabel(r'$M-M_{\odot}/h$')
        ax2.plot(r_bcg, m_s[k,:], 'r-', label=r'$M_\ast$')
        ax2.legend(loc = 2)
        ax2.set_xlim(1e-1, np.max(r_bcg))
        ax2.set_xticklabels(ax2.get_xticks(),rotation = 45,fontsize = 7.5)
        ax2.set_xscale('log')
        ax2.set_xlabel(r'$r-kpc/h$')
        ax2.set_title(r'$BCG m-r_{%.3f} MU_{%.3f}$' % (zmu[k], resolution))
        ax3.plot(r_bcg, m_b[k,:], 'b-', label=r'$M_b$')
        ax3.legend(loc = 2)
        ax3.set_xlim(1e-1, np.max(r_bcg))
        ax3.set_xticklabels(ax3.get_xticks(),rotation = 45,fontsize = 7.5)
        ax3.set_xscale('log')
        ax3.set_xlabel(r'$r-kpc/h$')
        plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
        plt.savefig(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_m-r_MU_No_%.0f_z_%.0f.png' % (_id_, k), dpi=600)
        plt.show()
        plt.close()
        '''
        plt.plot(r_bcg, m_g[k,:], 'g-', label=r'$M_g$')
        plt.plot(r_bcg, m_s[k,:], 'r-', label=r'$M_\ast$')
        plt.plot(r_bcg, m_b[k,:], 'b-', label=r'$M_b$')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e8,1e13)
        plt.xlim(1e-1, np.max(r_bcg))
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$M-M_{\odot}/h$')
        plt.legend(loc=2)
        plt.title(r'$BCG m-r_{%.3f} MU_{%.3f}$' % (zmu[k], resolution))
        plt.savefig(
        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_m-r_MU_No_%.0f_z_%.0f.png' % (_id_, k), dpi=600)
        plt.show()
        plt.close()
        '''
        plt.figure()
        f, (bx1, bx2, bx3) = plt.subplots(1,3,sharex = False,sharey=True)
        bx1.plot(r_bcg[:-1], rho_g[k,:], 'g-', label=r'$\rho_g$')
        bx1.legend(loc = 1)
        bx1.set_ylim(1e5,1e11)
        bx1.set_xlim(1e-1,np.max(r_bcg))
        bx1.set_xticklabels(bx1.get_xticks(),rotation = 45,fontsize = 7.5)
        bx1.set_yticklabels(bx1.get_yticks(),rotation = 90,fontsize = 7.5)
        bx1.set_xscale('log')
        bx1.set_yscale('log') 
        bx1.set_xlabel(r'$r-kpc/h$')
        bx1.set_ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')  
        bx2.plot(r_bcg[:-1], rho_s[k,:], 'r-', label=r'$\rho_{\ast}$')
        bx2.legend(loc = 1)
        bx2.set_xlim(1e-1,np.max(r_bcg))
        bx2.set_xticklabels(bx2.get_xticks(),rotation = 45,fontsize = 7.5)
        bx2.set_xscale('log')
        bx2.set_xlabel(r'$r-kpc/h$') 
        bx2.set_title(r'$BCG \rho-r_{%.3f} MU_{%.3f}$' % (zmu[k], resolution))
        bx3.plot(r_bcg[:-1], rho_b[k,:], 'b-', label=r'$\rho_b$')
        bx3.legend(loc = 1)
        bx3.set_xlim(1e-1,np.max(r_bcg))
        bx3.set_xticklabels(bx3.get_xticks(),rotation = 45,fontsize = 7.5)
        bx3.set_xscale('log')
        bx3.set_xlabel(r'$r-kpc/h$')
        plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)  
        plt.savefig(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_rho-r_MU_No_%.0f_z_%.0f.png' % (_id_, k), dpi=600)          
        plt.show()
        plt.close()
        '''
        plt.plot(r_bcg[:-1], rho_g[k,:], 'g-', label=r'$\rho_g$')
        plt.plot(r_bcg[:-1], rho_s[k,:], 'r-', label=r'$\rho_{\ast}$')
        plt.plot(r_bcg[:-1], rho_b[k,:], 'b-', label=r'$\rho_b$')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e5,1e11)
        plt.xlim(1e-1, np.max(r_bcg))
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')
        plt.legend(loc=1)
        plt.title(r'$BCG \rho-r_{%.3f} MU_{%.3f}$' % (zmu[k], resolution))
        plt.savefig(
        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_rho-r_MU_No_%.0f_z_%.0f.png' % (_id_, k), dpi=600)
        plt.show()
        plt.close()
        '''
        plt.figure()
        f, (cx1, cx2, cx3) = plt.subplots(1,3,sharex = False,sharey=True)
        cx1.plot(r_bcg, m_rho_g[k,:], 'g-', label=r'$\bar{\rho_g} MU$')
        cx1.legend(loc = 1)
        cx1.set_ylim(1e5,1e11)
        cx1.set_xlim(1e-1, np.max(r_bcg))
        cx1.set_xticklabels(cx1.get_xticks(),rotation = 45,fontsize = 7.5)
        cx1.set_yticklabels(cx1.get_yticks(),rotation = 90,fontsize = 7.5)
        cx1.set_xscale('log')
        cx1.set_yscale('log')
        cx1.set_xlabel(r'$r-kpc/h$')
        cx1.set_ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')  
        cx2.plot(r_bcg, m_rho_s[k,:], 'r-', label=r'$\bar{\rho_{\ast}} MU$')
        cx2.legend(loc = 1)
        cx2.set_xlim(1e-1, np.max(r_bcg))
        cx2.set_xticklabels(cx2.get_xticks(),rotation = 45,fontsize = 7.5)
        cx2.set_xscale('log')
        cx2.set_xlabel(r'$r-kpc/h$')
        cx2.set_title(r'$BCG \bar{\rho}-r_{%.3f} MU_{%.3f}$' % (zmu[k], resolution))
        cx3.plot(r_bcg, m_rho_b[k,:], 'b-', label=r'$\bar{\rho_b} MU$')
        cx3.legend(loc = 1)
        cx3.set_xlim(1e-1, np.max(r_bcg))
        cx3.set_xticklabels(cx3.get_xticks(),rotation = 45,fontsize = 7.5)
        cx3.set_xscale('log')
        cx3.set_xlabel(r'$r-kpc/h$')
        plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1) 
        plt.savefig(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_mean_rho-r_MU_No_%.0f_z_%.0f.png' % (_id_, k), dpi=600)            
        plt.show()
        plt.close() 
        '''
        plt.plot(r_bcg, m_rho_g[k,:], 'g-', label=r'$\bar{\rho_g} MU$')
        plt.plot(r_bcg, m_rho_s[k,:], 'r-', label=r'$\bar{\rho_{\ast}} MU$')
        plt.plot(r_bcg, m_rho_b[k,:], 'b-', label=r'$\bar{\rho_b} MU$')
        plt.legend(loc=1)
        plt.ylim(1e5,1e11)
        plt.xlim(1e-1, np.max(r_bcg))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')
        plt.title(r'$BCG \bar{\rho}-r_{%.3f} MU_{%.3f}$' % (zmu[k], resolution))
        plt.savefig(
        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_mean_rho-r_MU_No_%.0f_z_%.0f.png' % (_id_, k), dpi=600)
        plt.show()
        plt.close()
        '''
    for p in range(len(zgx)):
    #figure the GX simulation
        plt.figure()
        f, (ax1_gx, ax2_gx, ax3_gx) = plt.subplots(1,3,sharex = False,sharey=True)
        ax1_gx.plot(r_bcg_gx,m_g_gx[p,:],'g-',label = r'$M_g$')
        ax1_gx.legend(loc = 2)
        ax1_gx.set_ylim(1e8,1e13)
        ax1_gx.set_xlim(1e-1, np.max(r_bcg_gx))
        ax1_gx.set_xticklabels(ax1_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        ax1_gx.set_yticklabels(ax1_gx.get_yticks(),rotation = 90,fontsize = 7.5)
        ax1_gx.set_yscale('log')
        ax1_gx.set_xscale('log')
        ax1_gx.set_xlabel(r'$r-kpc/h$')
        ax1_gx.set_ylabel(r'$M-M_{\odot}/h$')
        ax2_gx.plot(r_bcg_gx,m_s_gx[p,:],'r-',label = r'$M_\ast$')
        ax2_gx.legend(loc = 2)
        ax2_gx.set_xlim(1e-1, np.max(r_bcg_gx))
        ax2_gx.set_xticklabels(ax2_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        ax2_gx.set_xscale('log')
        ax2_gx.set_xlabel(r'$r-kpc/h$')
        ax2_gx.set_title(r'$BCG m-r_{%.3f} GX_{%.3f}$'%(zgx[p],resolution))
        ax3_gx.plot(r_bcg_gx,m_b_gx[p,:],'b-',label = r'$M_b$')
        ax3_gx.legend(loc = 2)
        ax3_gx.set_xlim(1e-1, np.max(r_bcg_gx))
        ax3_gx.set_xticklabels(ax3_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        ax3_gx.set_xscale('log')
        ax3_gx.set_xlabel(r'$r-kpc/h$')
        plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
        plt.savefig(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_m-r_GX_No_%.0f_z_%.0f.png'%(_id_,p),dpi = 600)
        plt.show()
        plt.close()
        '''
        plt.plot(r_bcg_gx,m_g_gx[p,:],'g-',label = r'$M_g$')
        plt.plot(r_bcg_gx,m_s_gx[p,:],'r-',label = r'$M_\ast$')
        plt.plot(r_bcg_gx,m_b_gx[p,:],'b-',label = r'$M_b$')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e8,1e13)
        plt.xlim(1e-1,np.max(r_bcg_gx))
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$M-M_{\odot}/h$')
        plt.legend(loc = 2)
        plt.title(r'$BCG m-r_{%.3f} GX_{%.3f}$'%(zgx[p],resolution))
        plt.savefig(
        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_m-r_GX_No_%.0f_z_%.0f.png'%(_id_,p),dpi = 600)
        plt.show()
        plt.close()
        '''
        plt.figure()
        f, (bx1_gx, bx2_gx, bx3_gx) = plt.subplots(1,3,sharex = False,sharey=True)
        bx1_gx.plot(r_bcg_gx[:-1],rho_g_gx[p,:],'g-',label = r'$\rho_g$')
        bx1_gx.legend(loc = 1)
        bx1_gx.set_ylim(1e5,1e11)
        bx1_gx.set_xlim(1e-1,np.max(r_bcg_gx))
        bx1_gx.set_xticklabels(bx1_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        bx1_gx.set_yticklabels(bx1_gx.get_yticks(),rotation = 90,fontsize = 7.5)
        bx1_gx.set_xscale('log')
        bx1_gx.set_yscale('log') 
        bx1_gx.set_xlabel(r'$r-kpc/h$')
        bx1_gx.set_ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')  
        bx2_gx.plot(r_bcg_gx[:-1],rho_s_gx[p,:],'r-',label = r'$\rho_{\ast}$')
        bx2_gx.legend(loc = 1)
        bx2_gx.set_xlim(1e-1,np.max(r_bcg_gx))
        bx2_gx.set_xticklabels(bx2_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        bx2_gx.set_xscale('log')
        bx2_gx.set_xlabel(r'$r-kpc/h$') 
        bx2_gx.set_title(r'$BCG \rho-r_{%.3f} GX_{%.3f}$'%(zgx[p],resolution))
        bx3_gx.plot(r_bcg_gx[:-1],rho_b_gx[p,:],'b-',label = r'$\rho_b$')
        bx3_gx.legend(loc = 1)
        bx3_gx.set_xlim(1e-1,np.max(r_bcg_gx))
        bx3_gx.set_xticklabels(bx3_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        bx3_gx.set_xscale('log')
        bx3_gx.set_xlabel(r'$r-kpc/h$')
        plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)  
        plt.savefig(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_rho-r_GX_No_%.0f_z_%.0f.png'%(_id_,p),dpi = 600)
        plt.show()
        plt.close()
        '''
        plt.plot(r_bcg_gx[:-1],rho_g_gx[p,:],'g-',label = r'$\rho_g$')
        plt.plot(r_bcg_gx[:-1],rho_s_gx[p,:],'r-',label = r'$\rho_{\ast}$')
        plt.plot(r_bcg_gx[:-1],rho_b_gx[p,:],'b-',label = r'$\rho_b$')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e5,1e11)
        plt.xlim(1e-1,np.max(r_bcg_gx))
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')
        plt.legend(loc = 1)
        plt.title(r'$BCG \rho-r_{%.3f} GX_{%.3f}$'%(zgx[p],resolution))
        plt.savefig(
        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_rho-r_GX_No_%.0f_z_%.0f.png'%(_id_,p),dpi = 600)
        plt.show()
        plt.close()
        '''
        plt.figure()
        f, (cx1_gx, cx2_gx, cx3_gx) = plt.subplots(1,3,sharex = False,sharey=True)
        cx1_gx.plot(r_bcg_gx,m_rho_g_gx[p,:],'g-',label = r'$\bar{\rho_g}$')
        cx1_gx.legend(loc = 1)
        cx1_gx.set_ylim(1e5,1e11)
        cx1_gx.set_xlim(1e-1, np.max(r_bcg_gx))
        cx1_gx.set_xticklabels(cx1_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        cx1_gx.set_yticklabels(cx1_gx.get_yticks(),rotation = 90,fontsize = 7.5)
        cx1_gx.set_xscale('log')
        cx1_gx.set_yscale('log')
        cx1_gx.set_xlabel(r'$r-kpc/h$')
        cx1_gx.set_ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')  
        cx2_gx.plot(r_bcg_gx,m_rho_s_gx[p,:],'r-',label = r'$\bar{\rho_{\ast}}$')
        cx2_gx.legend(loc = 1)
        cx2_gx.set_xlim(1e-1, np.max(r_bcg_gx))
        cx2_gx.set_xticklabels(cx2_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        cx2_gx.set_xscale('log')
        cx2_gx.set_xlabel(r'$r-kpc/h$')
        cx2_gx.set_title(r'$BCG \bar{\rho}-r_{%.3f} GX_{%.3f}$'%(zgx[p],resolution))
        cx3_gx.plot(r_bcg_gx,m_rho_b_gx[p,:],'b-', label = r'$\bar{\rho_b}$')
        cx3_gx.legend(loc = 1)
        cx3_gx.set_xlim(1e-1, np.max(r_bcg_gx))
        cx3_gx.set_xticklabels(cx3_gx.get_xticks(),rotation = 45,fontsize = 7.5)
        cx3_gx.set_xscale('log')
        cx3_gx.set_xlabel(r'$r-kpc/h$')
        plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1) 
        plt.savefig(
                '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_mean_rho-r_GX_No_%.0f_z_%.0f.png'%(_id_,p),dpi = 600)           
        plt.show()
        plt.close()
        '''
        plt.plot(r_bcg_gx,m_rho_g_gx[p,:],'g-',label = r'$\bar{\rho_g}$')
        plt.plot(r_bcg_gx,m_rho_s_gx[p,:],'r-',label = r'$\bar{\rho_{\ast}}$')
        plt.plot(r_bcg_gx,m_rho_b_gx[p,:],'b-', label = r'$\bar{\rho_b}$')
        plt.legend(loc =1)
        plt.ylim(1e5,1e11)
        plt.xlim(1e-1,np.max(r_bcg))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')
        plt.title(r'$BCG \bar{\rho}-r_{%.3f} GX_{%.3f}$'%(zgx[p],resolution))
        plt.savefig(
        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_mean_rho-r_GX_No_%.0f_z_%.0f.png'%(_id_,p),dpi = 600)
        plt.show()
        plt.close()
        '''
    N = len(zmu)
    M = len(zgx)
    L_z = np.min([N,M])
    for q in range(L_z):
        #figure the comparation
        if zmu[q] <=id_z and zgx[q] <=id_z_gx:
            plt.figure()
            f,(kx1,kx2,kx3) = plt.subplots(1,3,sharex = False,sharey=True)
            kx1.plot(r_bcg,m_g[q,:],'g-',label = r'$M_g MU$')
            kx1.plot(r_bcg_gx,m_g_gx[q,:],'g--',label = r'$M_g GX$')
            kx1.legend(loc = 2)
            kx1.set_ylim(1e8,1e13)
            kx1.set_xlim(1e-1,np.max(r_bcg))
            kx1.set_xticklabels(kx1.get_xticks(),rotation = 45,fontsize = 7.5)
            kx1.set_yticklabels(kx1.get_yticks(),rotation = 90,fontsize = 7.5)
            kx1.set_xscale('log')
            kx1.set_yscale('log')
            kx1.set_xlabel(r'$r-kpc/h$')
            kx1.set_ylabel(r'$M-M_{\odot}/h$')
            kx2.plot(r_bcg,m_s[q,:],'r-',label = r'$M_\ast MU$')
            kx2.plot(r_bcg_gx,m_s_gx[q,:],'r--',label = r'$M_\ast GX$')
            kx2.legend(loc = 2)
            kx2.set_xlim(1e-1,np.max(r_bcg))
            kx2.set_xticklabels(kx2.get_xticks(),rotation = 45,fontsize = 7.5)
            kx2.set_xscale('log')
            kx2.set_xlabel(r'$r-kpc/h$')
            kx2.set_title(r'$BCG m-r_{%.3f} Re_{%.3f}$'%(zmu[q],resolution))
            kx3.plot(r_bcg,m_b[q,:],'b-',label = r'$M_b MU$')
            kx3.plot(r_bcg_gx,m_b_gx[q,:],'b--',label = r'$M_b GX$')
            kx3.legend(loc = 2)
            kx3.set_xlim(1e-1,np.max(r_bcg))
            kx3.set_xticklabels(kx3.get_xticks(),rotation = 45,fontsize = 7.5)
            kx3.set_xscale('log')
            kx3.set_xlabel(r'$r-kpc/h$')   
            plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1) 
            plt.savefig(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_m-r_No_%.0f_z_%.0f.png'%(_id_,q), dpi = 600)
            plt.show()
            plt.close()
            '''
            plt.plot(r_bcg,m_g[q,:],'g-',label = r'$M_g MU$')
            plt.plot(r_bcg_gx,m_g_gx[q,:],'g--',label = r'$M_g GX$')  
            plt.plot(r_bcg,m_s[q,:],'r-',label = r'$M_\ast MU$')
            plt.plot(r_bcg_gx,m_s_gx[q,:],'r--',label = r'$M_\ast GX$')
            plt.plot(r_bcg,m_b[q,:],'b-',label = r'$M_b MU$')
            plt.plot(r_bcg_gx,m_b_gx[q,:],'b--',label = r'$M_b GX$')
            plt.xscale('log')
            plt.yscale('log')
            plt.ylim(1e8,1e13)
            plt.xlim(1e-1,np.max(r_bcg))
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$M-M_{\odot}/h$')
            plt.legend(loc = 2)
            plt.title(r'$BCG m-r_{%.3f} Re_{%.3f}$'%(zmu[q],resolution))
            plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_m-r_No_%.0f_z_%.0f.png'%(_id_,q), dpi = 600)
            plt.show()
            plt.close()
            '''
            plt.figure()
            f,(lx1,lx2,lx3) = plt.subplots(1,3,sharex = False,sharey=True)
            lx1.plot(r_bcg[:-1],rho_g[q,:],'g-',label = r'$\rho_g MU$')
            lx1.plot(r_bcg_gx[:-1],rho_g_gx[q,:],'g--',label = r'$\rho_g GX$')
            lx1.legend(loc = 1)
            lx1.set_ylim(1e5,1e11)
            lx1.set_xlim(1e-1,np.max(r_bcg))
            lx1.set_xticklabels(lx1.get_xticks(),rotation = 45,fontsize = 7.5)
            lx1.set_yticklabels(lx1.get_yticks(),rotation = 90,fontsize = 7.5)
            lx1.set_xscale('log')
            lx1.set_yscale('log')
            lx1.set_xlabel(r'$r-kpc/h$')
            lx1.set_ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')  
            lx2.plot(r_bcg[:-1],rho_s[q,:],'r-',label = r'$\rho_{\ast} MU$')
            lx2.plot(r_bcg_gx[:-1],rho_s_gx[q,:],'r--',label = r'$\rho_{\ast} GX$')
            lx2.legend(loc = 1)
            lx2.set_xlim(1e-1,np.max(r_bcg))
            lx2.set_xticklabels(lx2.get_xticks(),rotation = 45,fontsize = 7.5)
            lx2.set_xscale('log') 
            lx2.set_xlabel(r'$r-kpc/h$')  
            lx2.set_title(r'$BCG \rho-r_{%.3f} Re_{%.3f}$'%(zmu[q],resolution))
            lx3.plot(r_bcg[:-1],rho_b[q,:],'b-',label = r'$\rho_b MU$')
            lx3.plot(r_bcg_gx[:-1],rho_b_gx[q,:],'b--',label = r'$\rho_b GX$')
            lx3.legend(loc = 1)
            lx3.set_xlim(1e-1,np.max(r_bcg))
            lx3.set_xticklabels(lx3.get_xticks(),rotation = 45,fontsize = 7.5)
            lx3.set_xscale('log') 
            lx3.set_xlabel(r'$r-kpc/h$')  
            plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1) 
            plt.savefig(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_rho-r_No_%.0f_z_%.0f.png'%(_id_,q), dpi = 600)
            plt.show()
            plt.close()
            '''
            plt.plot(r_bcg[:-1],rho_g[q,:],'g-',label = r'$\rho_g MU$')
            plt.plot(r_bcg[:-1],rho_s[q,:],'r-',label = r'$\rho_{\ast} MU$')
            plt.plot(r_bcg[:-1],rho_b[q,:],'b-',label = r'$\rho_b MU$')
            plt.plot(r_bcg_gx[:-1],rho_g_gx[q,:],'g--',label = r'$\rho_g GX$')
            plt.plot(r_bcg_gx[:-1],rho_s_gx[q,:],'r--',label = r'$\rho_{\ast} GX$')
            plt.plot(r_bcg_gx[:-1],rho_b_gx[q,:],'b--',label = r'$\rho_b GX$')
            plt.xscale('log')
            plt.yscale('log')
            plt.ylim(1e5,1e11)
            plt.xlim(1e-1,np.max(r_bcg))
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')
            plt.legend(loc = 1)
            plt.title(r'$BCG \rho-r_{%.3f} Re_{%.3f}$'%(zmu[q],resolution))
            plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_rho-r_No_%.0f_z_%.0f.png'%(_id_,q), dpi = 600)
            plt.show()
            plt.close()
            '''
            plt.figure()
            f,(mx1,mx2,mx3) = plt.subplots(1,3,sharex = False,sharey=True)
            mx1.plot(r_bcg,m_rho_g[q,:],'g-',label = r'$\bar{\rho_g} MU$')
            mx1.plot(r_bcg_gx,m_rho_g_gx[q,:],'g--',label = r'$\bar{\rho_g} GX$')
            mx1.legend(loc = 1)
            mx1.set_ylim(1e5,1e11)
            mx1.set_xlim(1e-1,np.max(r_bcg))
            mx1.set_xticklabels(mx1.get_xticks(),rotation = 45,fontsize = 7.5)
            mx1.set_yticklabels(mx1.get_yticks(),rotation = 90,fontsize = 7.5)
            mx1.set_xscale('log')
            mx1.set_yscale('log')
            mx1.set_xlabel(r'$r-kpc/h$')
            mx1.set_ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')
            mx2.plot(r_bcg,m_rho_s[q,:],'r-',label = r'$\bar{\rho_{\ast}} MU$')
            mx2.plot(r_bcg_gx,m_rho_s_gx[q,:],'r--',label = r'$\bar{\rho_{\ast}} GX$')
            mx2.legend(loc = 1)
            mx2.set_xlim(1e-1,np.max(r_bcg))
            mx2.set_xticklabels(mx2.get_xticks(),rotation = 45,fontsize = 7.5)
            mx2.set_xscale('log')            
            mx2.set_xlabel(r'$r-kpc/h$')  
            mx2.set_title(r'$BCG \bar{\rho}-r_{%.3f} Re_{%.3f}$'%(zmu[q],resolution))
            mx3.plot(r_bcg,m_rho_b[q,:],'b-',label = r'$\bar{\rho_b} MU$')
            mx3.plot(r_bcg_gx,m_rho_b_gx[q,:],'b--',label = r'$\bar{\rho_b} GX$')
            mx3.legend(loc = 1)
            mx3.set_xlim(1e-1,np.max(r_bcg))
            mx3.set_xticklabels(mx3.get_xticks(),rotation = 45,fontsize = 7.5)
            mx3.set_xscale('log')            
            mx3.set_xlabel(r'$r-kpc/h$')
            plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1) 
            plt.savefig(
                    '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_mean-rho_No_%.0f_z_%.0f.png'%(_id_,q), dpi = 600)
            plt.show()
            plt.close()
            '''
            plt.plot(r_bcg,mean_rho_g,'g-',label = r'$\bar{\rho_g} MU$')
            plt.plot(r_bcg,mean_rho_s,'r-',label = r'$\bar{\rho_{\ast}} MU$')
            plt.plot(r_bcg,mean_rho_b,'b-',label = r'$\bar{\rho_b} MU$')
            plt.plot(r_bcg_gx,mean_rho_g_gx,'g--',label = r'$\bar{\rho_g} GX$')
            plt.plot(r_bcg_gx,mean_rho_s_gx,'r--',label = r'$\bar{\rho_{\ast}} GX$')
            plt.plot(r_bcg_gx,mean_rho_b_gx,'b--',label = r'$\bar{\rho_b} GX$')
            plt.legend(loc = 1)
            plt.ylim(1e5,1e11)
            plt.xlim(1e-1,np.max(r_bcg))
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')
            plt.title(r'$BCG \bar{\rho}-r_{%.3f} Re_{%.3f}$'%(zmu[q],resolution))
            plt.savefig(
            '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/snap/BCG_fig/BCG_mean-rho_No_%.0f_z_%.0f.png'%(_id_,q), dpi = 600)
            plt.show()
            plt.close()
            '''
    return
#find_BCG_fig(ip = True,resolution = True,g_size = True)