"""
###this file try to figure out all of the mass evolution of main halo
##from the redshift of the merger start
"""
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
import changds
#import glob
#def tot_halo_profile(ip):
_id_ = 0
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)
#read the meger tree to get the star time of each halo
with h5py.File('/home/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/Redshift.h5') as f:
    y0 = f['a']
    com_z = np.array(y0)
with h5py.File('/home/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/main_tree.h5') as f:
    y1 = f['a']
    tree_line = np.array(y1)
with h5py.File('/home/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/Redshift_GX.h5') as f:
    y2 = f['a']
    com_z_gx = np.array(y2)
with h5py.File('/home/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/main_tree_GX.h5') as f:
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
###set array use to figure
rhog = {}
rhos = {}
rhob = {}
rhod = {}
mrhog = {}
mrhos = {}
mrhob = {}
mrhod = {}
R = {}
zmu = {}
rhog_gx = {}
rhos_gx = {}
rhob_gx = {}
rhod_gx = {}
mrhog_gx = {}
mrhos_gx = {}
mrhob_gx = {}
mrhod_gx = {}
R_gx = {}
zgx = {}
for k in range(len(alpha)):
    dd = alpha[-k]
    ss = str(dd)
    if len(ss) == 1:
        s_u = str('00%.0f' % dd)
    elif len(ss) == 2:
        s_u = str('0%.0f' % dd)
    else:
        s_u = str(dd)
    No_snap = s_u
    snap_z = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
    snap_z = np.abs(snap_z)
    print('Now redshift is %.3f' % snap_z) 
    if snap_z <= id_z:
        zmu[k] = snap_z
        goal_z = changds.inv_chidas('%.3f'%snap_z)
        try:
            with h5py.File('Read_halo-%.0f-profile_MU-%.3f.h5'%(k,goal_z)) as f:
                m_MU = np.array(f['a'])
            mg = m_MU[0,:]
            ms = m_MU[1,:]
            mb = m_MU[2,:]
            md = m_MU[3,:]
            r = m_MU[4,:]
            M = len(r)
            rho_g = np.zeros(M-1,dtype = np.float)
            rho_s = np.zeros(M-1,dtype = np.float)
            rho_b = np.zeros(M-1,dtype = np.float)
            rho_d = np.zeros(M-1,dtype = np.float)
            for t in range(M-1):
                dr = r[t+1] - r[t]
                rho_g[t] = (mg[t+1]-mg[t])/(4*np.pi*r[t]**2*dr)
                rho_s[t] = (ms[t+1]-ms[t])/(4*np.pi*r[t]**2*dr)
                rho_b[t] = (mb[t+1]-mb[t])/(4*np.pi*r[t]**2*dr)
                rho_d[t] = (md[t+1]-md[t])/(4*np.pi*r[t]**2*dr)
            rhog[k] = rho_g
            rhos[k] = rho_s
            rhob[k] = rho_b
            rhod[k] = rho_d
            R[k] = r
            mrhog[k] = mg/(4*np.pi*r**3/3)
            mrhos[k] = ms/(4*np.pi*r**3/3)
            mrhob[k] = mb/(4*np.pi*r**3/3)
            mrhod[k] = md/(4*np.pi*r**3/3)
        except KeyError:
            print('Now have read out all of the data!')
            break
for k in range(len(alpha)):
    dd = alpha[-k]
    ss = str(dd)
    if len(ss) == 1:
        s_u = str('00%.0f' % dd)
    elif len(ss) == 2:
        s_u = str('0%.0f' % dd)
    else:
        s_u = str(dd)
    No_snap = s_u
    snap_z_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
    snap_z_gx = np.abs(snap_z_gx)
    print('Now redshift is %.3f' % snap_z_gx)
    if snap_z_gx <= id_z_gx:
        zgx[k] = snap_z_gx
        goal_z_gx = changds.inv_chidas('%.3f'%snap_z_gx)
        try:
            with h5py.File('Read_halo-%.0f-profile_GX-%.3f.h5'%(k,goal_z_gx)) as f:
                m_GX = np.array(f['a'])
            mg_gx = m_GX[0,:]
            ms_gx = m_GX[1,:]
            mb_gx = m_GX[2,:]
            md_gx = m_GX[3,:]
            r_gx = m_GX[4,:]
            N = len(r_gx)
            rho_g_gx = np.zeros(N-1,dtype = np.float)
            rho_s_gx = np.zeros(N-1,dtype = np.float)
            rho_b_gx = np.zeros(N-1,dtype = np.float)
            rho_d_gx = np.zeros(N-1,dtype = np.float)
            for t in range(N-1):
                dr_gx = r_gx[t+1] - r_gx[t]
                rho_g_gx[t] = (mg_gx[t+1]-mg_gx[t])/(4*np.pi*r_gx[t]**2*dr_gx)
                rho_s_gx[t] = (ms_gx[t+1]-ms_gx[t])/(4*np.pi*r_gx[t]**2*dr_gx)
                rho_b_gx[t] = (mb_gx[t+1]-mb_gx[t])/(4*np.pi*r_gx[t]**2*dr_gx)
                rho_d_gx[t] = (md_gx[t+1]-md_gx[t])/(4*np.pi*r_gx[t]**2*dr_gx)
            rhog_gx[k] = rho_g_gx
            rhos_gx[k] = rho_s_gx
            rhob_gx[k] = rho_b_gx
            rhod_gx[k] = rho_d_gx        
            R_gx[k] = r_gx
            mrhog_gx[k] = mg_gx/(4*np.pi*r_gx**3/3)
            mrhos_gx[k] = ms_gx/(4*np.pi*r_gx**3/3)
            mrhob_gx[k] = mb_gx/(4*np.pi*r_gx**3/3)
            mrhod_gx[k] = md_gx/(4*np.pi*r_gx**3/3)
        except KeyError:
            print('Now have read out all of the data!')
            break            
#figure out the results of MUSIC
f,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,sharex = False,sharey = True,figsize = (8,8))
L = len(R)
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        ax1.plot(r[:-1],rhog[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
ax1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
ax1.set_xlim(1e0,1e3)
ax1.set_ylim(1e4,1e11)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
ax1.set_xticklabels(ax1.get_xticks(),rotation = 45,fontsize = 7.5)
ax1.set_yticklabels(ax1.get_yticks(),rotation = 90,fontsize = 7.5) 
plt.rcParams['figure.figsize'] = (5,5)
ax1.set_xlabel(r'$r [kpc/h]$')
ax1.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_title(r'$Halo-\rho_g-r-MU_{%.0f}$'%_id_,fontsize = 8.5)
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        ax2.plot(r[:-1],rhos[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
ax2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
ax2.set_xlim(1e0,1e3)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
ax2.set_xticklabels(ax2.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['figure.figsize'] = (5,5)
ax2.set_xscale('log')
ax2.set_xlabel(r'$r [kpc/h]$')
#ax2.set_ylim(1e4,1e11)
#ax2.set_yticklabels(ax2.get_yticks(),rotation = 90,fontsize = 7.5) 
#ax2.set_yscale('log')
#ax2.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
ax2.set_title(r'$Halo-\rho_{\ast}-r-MU_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        ax3.plot(r[:-1],rhob[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
ax3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
ax3.set_xlim(1e0,1e3)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in'
ax3.set_xticklabels(ax3.get_xticks(),rotation = 45,fontsize = 7.5) 
plt.rcParams['figure.figsize'] = (5,5)
ax3.set_xscale('log')
ax3.set_xlabel(r'$r [kpc/h]$')
#ax3.set_ylim(1e4,1e11)
#ax3.set_yticklabels(ax3.get_yticks(),rotation = 90,fontsize = 7.5) 
#ax3.set_yscale('log')
#ax3.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
ax3.set_title(r'$Halo-\rho_b-r-MU_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        ax4.plot(r[:-1],rhod[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
ax4.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
ax4.set_xlim(1e0,1e3)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in'
ax4.set_xticklabels(ax4.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['figure.figsize'] = (5,5)
ax4.set_xlabel(r'$r [kpc/h]$')
ax4.set_xscale('log')
#ax4.set_ylim(1e4,1e11)
#ax4.set_yticklabels(ax4.get_yticks(),rotation = 90,fontsize = 7.5) 
#ax4.set_yscale('log')
#ax4.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
ax4.set_title(r'$Halo-\rho_d-r-MU_{%.0f}$'%_id_,fontsize = 8.5) 
plt.tight_layout(pad=0.5, h_pad=0., w_pad= -1)
plt.savefig('Halo density profile evolution MU',dpi = 600)
plt.show()
plt.close()
#####
f,(cx1,cx2,cx3,cx4) = plt.subplots(1,4,sharex = False,sharey = True,figsize = (8,8))
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        cx1.plot(r,mrhog[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
cx1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
cx1.set_xlabel(r'$r [kpc/h]$')
cx1.set_xlim(1e0,1e3)
cx1.set_ylim(1e4,1e11)
cx1.set_xticklabels(cx1.get_xticks(),rotation = 45,fontsize = 7.5)
cx1.set_yticklabels(cx1.get_yticks(),rotation = 90,fontsize = 7.5) 
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
cx1.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
cx1.set_yscale('log')
cx1.set_xscale('log')
cx1.set_title(r'$Halo-\bar{\rho_g}-r-MU_{%.0f}$'%_id_,fontsize = 8.5)
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        cx2.plot(r,mrhos[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
cx2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
cx2.set_xlabel(r'$r [kpc/h]$')
cx2.set_xlim(1e0,1e3)
cx2.set_xticklabels(cx2.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['figure.figsize'] = (5,5) 
cx2.set_xscale('log')
#cx2.set_ylim(1e4,1e11)
#cx2.set_yticklabels(cx2.get_yticks(),rotation = 90,fontsize = 7.5) 
#cx2.set_yscale('log')
#cx2.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
cx2.set_title(r'$Halo-\bar{\rho_{\ast}}-r-MU_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        cx3.plot(r,mrhob[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
cx3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
cx3.set_xlabel(r'$r [kpc/h]$')
cx3.set_xlim(1e0,1e3)
cx3.set_xticklabels(cx3.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
cx3.set_xscale('log')
#cx3.set_ylim(1e4,1e11)
#cx3.set_yticklabels(cx3.get_yticks(),rotation = 90,fontsize = 7.5) 
#cx3.set_yscale('log')
#cx3.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
cx3.set_title(r'$Halo-\bar{\rho_b}-r-MU_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,L):
    r = R[p]
    if p % 8 == 1:
        cx4.plot(r,mrhod[p],'-',color = mpl.cm.rainbow(p/L),
                 label = r'$z{%.3f}$'%zmu[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
cx4.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
cx4.set_xlabel(r'$r [kpc/h]$')
cx4.set_xlim(1e0,1e3)
cx4.set_xticklabels(cx4.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
cx4.set_xscale('log')
#cx4.set_ylim(1e4,1e11)
#cx4.set_yticklabels(cx4.get_yticks(),rotation = 90,fontsize = 7.5) 
#cx4.set_yscale('log')
#cx4.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
cx4.set_title(r'$Halo-\bar{\rho_d}-r-MU_{%.0f}$'%_id_,fontsize = 8.5) 
plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
plt.savefig('Halo mean density profile evolution MU',dpi = 600)
plt.show()
plt.close()
#figure out the results of GadgetX 
f,(bx1,bx2,bx3,bx4) = plt.subplots(1,4,sharex = False,sharey = True,figsize = (8,8))
S = len(R_gx) 
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        bx1.plot(r_gx[:-1],rhog_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
bx1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
bx1.set_xlabel(r'$r [kpc/h]$')
bx1.set_xlim(1e0,1e3)
bx1.set_ylim(1e4,1e11)
bx1.set_xticklabels(bx1.get_xticks(),rotation = 45,fontsize = 7.5)
bx1.set_yticklabels(bx1.get_yticks(),rotation = 90,fontsize = 7.5) 
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
bx1.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
bx1.set_yscale('log')
bx1.set_xscale('log')
bx1.set_title(r'$Halo-\rho_g-r-GX_{%.0f}$'%_id_,fontsize = 8.5)
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        bx2.plot(r_gx[:-1],rhos_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
bx2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
bx2.set_xlabel(r'$r [kpc/h]$')
bx2.set_xlim(1e0,1e3)
bx2.set_xticklabels(bx2.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['figure.figsize'] = (5,5) 
bx2.set_xscale('log')
#bx2.set_ylim(1e4,1e11)
#bx2.set_yticklabels(bx2.get_yticks(),rotation = 90,fontsize = 7.5) 
#bx2.set_yscale('log')
#bx2.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
bx2.set_title(r'$Halo-\rho_{\ast}-r-GX_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        bx3.plot(r_gx[:-1],rhob_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
bx3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
bx3.set_xlabel(r'$r [kpc/h]$')
bx3.set_xlim(1e0,1e3)
bx3.set_xticklabels(bx3.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
bx3.set_xscale('log')
#bx3.set_ylim(1e4,1e11)
#bx3.set_yticklabels(bx3.get_yticks(),rotation = 90,fontsize = 7.5) 
#bx3.set_yscale('log')
#bx3.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
bx3.set_title(r'$Halo-\rho_b-r-GX_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        bx4.plot(r_gx[:-1],rhod_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
bx4.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
bx4.set_xlabel(r'$r [kpc/h]$')
bx4.set_xlim(1e0,1e3)
bx4.set_xticklabels(bx4.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
bx4.set_xscale('log')
#bx4.set_ylim(1e4,1e11)
#bx4.set_yticklabels(bx4.get_yticks(),rotation = 90,fontsize = 7.5) 
#bx4.set_yscale('log')
#bx4.set_ylabel(r'$\rho[M_\odot h^2/kpc^3]$')
bx4.set_title(r'$Halo-\rho_d-r-MU_{%.0f}$'%_id_,fontsize = 8.5) 
plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
plt.savefig('Halo density profile evolution GX',dpi = 600)
plt.show()
plt.close() 

f,(dx1,dx2,dx3,dx4) = plt.subplots(1,4,sharex = False,sharey = True,figsize = (8,8)) 
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        dx1.plot(r_gx,mrhog_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
dx1.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
dx1.set_xlabel(r'$r [kpc/h]$')
dx1.set_xlim(1e0,1e3)
dx1.set_ylim(1e4,1e11)
dx1.set_xticklabels(dx1.get_xticks(),rotation = 45,fontsize = 7.5)
dx1.set_yticklabels(dx1.get_yticks(),rotation = 90,fontsize = 7.5) 
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
dx1.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
dx1.set_yscale('log')
dx1.set_xscale('log')
dx1.set_title(r'$Halo-\bar{\rho_g}-r-GX_{%.0f}$'%_id_,fontsize = 8.5)
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        dx2.plot(r_gx,mrhos_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
dx2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
dx2.set_xlabel(r'$r [kpc/h]$')
dx2.set_xlim(1e0,1e3)
dx2.set_xticklabels(dx2.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
dx2.set_xscale('log')
#dx2.set_ylim(1e4,1e11)
#dx2.set_yticklabels(dx2.get_yticks(),rotation = 90,fontsize = 7.5) 
#dx2.set_yscale('log')
#dx2.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
dx2.set_title(r'$Halo-\bar{\rho_{\ast}}-r-GX_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        dx3.plot(r_gx,mrhob_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
dx3.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
dx3.set_xlabel(r'$r [kpc/h]$')
dx3.set_xlim(1e0,1e3)
dx3.set_xticklabels(dx3.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
dx3.set_xscale('log')
#dx3.set_ylim(1e4,1e11)
#dx3.set_yticklabels(dx3.get_yticks(),rotation = 90,fontsize = 7.5) 
#dx3.set_yscale('log')
#dx3.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
dx3.set_title(r'$Halo-\bar{\rho_b}-r-GX_{%.0f}$'%_id_,fontsize = 8.5)   
for p in range(1,S):
    r_gx = R_gx[p]
    if p % 8 == 1:
        dx4.plot(r_gx,mrhod_gx[p],'-',color = mpl.cm.rainbow(p/S),
                 label = r'$z{%.3f}$'%zgx[p])
        handles, labels = plt.gca().get_legend_handles_labels() 
dx4.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol= 3, mode="expand", fontsize= 6.5, borderaxespad=0.)
dx4.set_xlabel(r'$r [kpc/h]$')
dx4.set_xlim(1e0,1e3)
dx4.set_xticklabels(dx4.get_xticks(),rotation = 45,fontsize = 7.5)
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in' 
plt.rcParams['figure.figsize'] = (5,5)
dx4.set_xscale('log')
#dx4.set_ylim(1e4,1e11)
#dx4.set_yticklabels(dx4.get_yticks(),rotation = 90,fontsize = 7.5) 
#dx4.set_yscale('log')
#dx4.set_ylabel(r'$\bar{\rho}[M_\odot h^2/kpc^3]$')
dx4.set_title(r'$Halo-\bar{\rho_d}-r-MU_{%.0f}$'%_id_,fontsize = 8.5) 
plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)
plt.savefig('Halo mean density profile evolution GX',dpi = 600)
plt.show()
plt.close() 
