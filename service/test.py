import matplotlib as mpl
mpl.use('Agg')
from pygadgetreader import *
import h5py
import scipy.stats as st
import find
import pygadgetreader as pygdr
import numpy as np
import astropy.io.ascii as asc
import matplotlib.pyplot as plt
import pandas as pd
import astropy.constants as C
import astropy.units as U
import itertools
import handy
import find 
#initial data
resolution = 1
size_BCG = 100
_id_ = 0
Nr = 1./resolution
# to make sure the bins density is the same
N_size = 50.*(size_BCG/100.0)
Nsize = np.int0(np.ceil(N_size))
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)
#read the meger tree data
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
ia1 = com_z <= id_z
z_mu = np.array(com_z[ia1])
ia2 = com_z_gx <= id_z_gx
z_gx = np.array(com_z_gx[ia2]) 
f, (ax1, ax2) = plt.subplots(1,2,sharex = False,sharey=True)
n_color = len(z_mu)
plot_lines1 = []
z_label1 = []
for k in range(len(z_mu)):
    try:
        if k %8 == 0:
            z_label1.append(z_mu[k])
            with h5py.File('BCG_T_R_h%.0f_z%.3f_MU.h5'%(_id_,z_mu[k])) as f:
                array5 = np.array(f['a'])
            r_bcg = np.array(array5[0,:])
            T_galaxy_m = np.array(array5[1,:])
            T_galaxy_n = np.array(array5[2,:])             
            l1, = ax1.plot(r_bcg,T_galaxy_m,'-',
                           color = mpl.cm.rainbow(k/n_color),label = r'$\bar{T}-r-in-M_{%.3f}$'%z_mu[k])
            l2, = ax1.plot(r_bcg,T_galaxy_n,'--',
                           color = mpl.cm.rainbow(k/n_color),label = r'$\bar{T}-r-in-N_{%.3f}$'%z_mu[k])
            plot_lines1.append([l1, l2])
            handles, labels = plt.gca().get_legend_handles_labels()
    except KeyError:
        break
#ax1.legend(bbox_to_anchor=(-0.05, 1., 1.15, 0.), loc=2,
#           ncol= 3, mode="expand", fontsize= 5.5, borderaxespad=2.)
#legend1 = plt.legend(plot_lines1[0],['average in mass','average in number'],loc = 2)
#plt.legend([l[0] for l in plot_lines1], z_label1, loc = 4)
#plt.gca().add_artist(legend1)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['figure.figsize'] = (5,5)
ax1.set_xlabel(r'$r-kpc/h$')
ax1.set_ylabel(r'$\bar{T}-K$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title(r'$T-R_{3D}-galaxy%.0f-MU$'%_id_,fontsize = 7.5)   

plot_lines2 = []
z_label2 = []
for k in range(len(z_mu)):
    try:
        if k %8 == 0:
            z_label2.append(z_mu[k])
            with h5py.File('meanT_R_h%.0f_z%.3f_MU.h5'%(_id_,z_mu[k])) as f:
                array4 = np.array(f['a'])
            R = np.array(array4[0,:])
            mean_T_m = np.array(array4[1,:])
            mean_T_n = np.array(array4[2,:])     
            l3, = ax2.plot(R,mean_T_m,'-',color = mpl.cm.rainbow(k/n_color),label = r'$\bar{T}-r-in-M_{%.3f}$'%z_mu[k])
            l4, = ax2.plot(R,mean_T_n,'--',color = mpl.cm.rainbow(k/n_color),label = r'$\bar{T}-r-in-N_{%.3f}$'%z_mu[k])   
            handles, labels = plt.gca().get_legend_handles_labels()
            plot_lines2.append([l3, l4])
    except KeyError:
        break
#ax2.legend(bbox_to_anchor=(-0.05, 1., 1.15, 0.), loc=2,
#           ncol= 3, mode="expand", fontsize= 5.5, borderaxespad=2.)
legend2 = plt.legend(plot_lines2[0],['average in mass','average in number'],loc = 2, fontsize= 7.5)
plt.legend([l[0] for l in plot_lines2], z_label2, loc = 4, fontsize= 7.5)
plt.gca().add_artist(legend2)
plt.rcParams['xtick.direction'] ='in'
plt.rcParams['ytick.direction'] ='in'
plt.rcParams['figure.figsize'] = (5,5)
ax2.set_xlabel(r'$r-kpc/h$')
#ax2.set_ylabel(r'$\bar{T}-K$')
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.set_title(r'$T-R_{3D}-halo%.0f-MU$'%_id_,fontsize = 7.5)
plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)     
plt.savefig('T_SFR/temperatur profile BCG MU%.0f'%_id_,dpi = 600)
plt.show()
plt.close()

f, (bx1, bx2) = plt.subplots(1,2,sharex = False,sharey=True)
m_color = len(z_gx)
plot_lines3 = []
z_label3 = []
for k in range(len(z_gx)):
    try:
        if k %8 == 0:
            z_label3.append(z_gx[k])
            with h5py.File('BCG_T_R_h%.0f_z%.3f_GX.h5'%(_id_, z_gx[k])) as f:
                array5_gx = np.array(f['a'])
            r_bcg_gx = np.array(array5_gx[0,:])
            T_galaxy_m_gx = np.array(array5_gx[1,:])
            T_galaxy_n_gx = np.array(array5_gx[2,:]) 
            l5, = bx1.plot(r_bcg_gx,T_galaxy_m_gx,'-',color = mpl.cm.rainbow(k/m_color),label = r'$\bar{T}-r-in-M_{%.3f}$'%z_gx[k])
            l6, = bx1.plot(r_bcg_gx,T_galaxy_n_gx,'--',color = mpl.cm.rainbow(k/m_color),label = r'$\bar{T}-r-in-N_{%.3f}$'%z_gx[k])
            handles, labels = plt.gca().get_legend_handles_labels()
            plot_lines3.append([l5, l6])
    except KeyError:
        break
#bx1.legend(bbox_to_anchor=(-0.05, 1., 1.15, 0.), loc=2,
#           ncol= 3, mode="expand", fontsize= 5.5, borderaxespad=2.)
#legend3 = plt.legend(plot_lines3[0],['average in mass','average in number'],loc = 2, fontsize= 7.5)
#plt.legend([l[0] for l in plot_lines3], z_label3, loc = 4, fontsize= 7.5)
#plt.gca().add_artist(legend3)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['figure.figsize'] = (5,5)
bx1.set_xlabel(r'$r-kpc/h$')
bx1.set_ylabel(r'$\bar{T}-K$')
bx1.set_xscale('log')
bx1.set_yscale('log')
bx1.set_title(r'$T-R_{3D}-galaxy%.0f-GX$'%_id_,fontsize = 7.5)

plot_lines4 = []
z_label4 = []
for k in range(len(z_gx)):
    try:
        if k %8 == 0:
            z_label4.append(z_gx[k])
            with h5py.File('meanT_R_h%.0f_z%.3f_GX.h5'%(_id_,z_gx[k])) as f:
                array4_gx = np.array(f['a'])    
            R_gx = np.array(array4_gx[0,:])
            mean_T_m_gx = np.array(array4_gx[1,:])
            mean_T_n_gx = np.array(array4_gx[2,:])       
            l7, = bx2.plot(R_gx,mean_T_m_gx,'-',color = mpl.cm.rainbow(k/n_color),label = r'$\bar{T}-r-in-M_{%.3f}$'%z_gx[k])
            l8, = bx2.plot(R_gx,mean_T_n_gx,'--',color = mpl.cm.rainbow(k/n_color),label = r'$\bar{T}-r-in-N_{%.3f}$'%z_gx[k])   
            handles, labels = plt.gca().get_legend_handles_labels()
            plot_lines4.append([l7, l8])
    except KeyError:
        break
#bx2.legend(bbox_to_anchor=(-0.05, 1., 1.15, 0.), loc=2,
#           ncol= 3, mode="expand", fontsize= 5.5, borderaxespad=2.)
legend4 = plt.legend(plot_lines4[0],['average in mass','average in number'],loc = 2, fontsize= 7.5)
plt.legend([l[0] for l in plot_lines4], z_label4, loc = 4, fontsize= 7.5)
plt.gca().add_artist(legend4)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['figure.figsize'] = (5,5)
bx2.set_xlabel(r'$r-kpc/h$')
#bx2.set_ylabel(r'$\bar{T}-K$')
bx2.set_xscale('log')
#bx2.set_yscale('log')
bx2.set_title(r'$T-R_{3D}-halo%.0f-GX$'%_id_,fontsize = 7.5)   
plt.tight_layout(pad=0.5, h_pad=0., w_pad=-1)     
plt.savefig('T_SFR/temperatur profile BCG GX%.0f'%_id_,dpi = 600)
plt.show()
plt.close()
