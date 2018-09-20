#this file try to find the BCG in a given halo(to find the most dengsity part of given halo)
##PS :the devided cut just for the central area of the given halo
from pygadgetreader import *
import pygadgetreader as pygdr
import numpy as np
import h5py
import pandas as pd
import astropy.io.ascii as asc
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib as mpl
from handy import scatter as hsc
import find
No_snap = '128'
_id_ = 0
resolution = 1.0
Nr = 1./resolution
size_BCG = 100
#to make sure the bins density is the same
N_size = 50.*(size_BCG/100.0)
N_size = np.ceil(N_size)
#set the goal snap(redshift),goal halo id,resolution(in nuit : kpc/h)
#set 1kpc/h as a measurement 
#######################section1:read out the position data of MUSIC
snap_z = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'redshift')
snap_z = np.abs(snap_z)
print('Now redshift is %.3f'%snap_z)
snap_N = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'npartTotal')
print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
      %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
snap_name = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'header')
#print(snap_name)
#for type
snap_shot_gas = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','gas')
snap_shot_DM = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','dm')
try:
    snap_shot_star = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','star')
    snap_mass_star = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'mass','star')
except SystemExit:
    print('no star particles now')
#for respective position
snap_shot_bulge = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','bulge')
snap_shot_disk = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','disk')
try:
    snap_shot_bndry = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','bndry')
    snap_mass_BH = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'mass','bndry')
except SystemExit:
    print('no boundary particles now')
snap_mass_gas = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'mass','gas')
#for density profile of gas(there only density about gas in the simulation data)
#snap_dens = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'rho','gas') 
#the total mass distribution of snap_shot_128
main_halo = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
Rvir = np.array(main_halo['col12'])
xhalo = np.array(main_halo['col6'])
yhalo = np.array(main_halo['col7'])
zhalo = np.array(main_halo['col8'])
x0 = xhalo[_id_]
y0 = yhalo[_id_]
z0 = zhalo[_id_]
R0 = Rvir[_id_]
dgas = np.sqrt((snap_shot_gas[:,0]-x0)**2+(snap_shot_gas[:,1]-y0)**2+
               (snap_shot_gas[:,2]-z0)**2)
ig = dgas <= R0
inlgas = snap_shot_gas[ig,:]
inlmass_gas = snap_mass_gas[ig]

if snap_N[-2] !=0:
   dstar = np.sqrt((snap_shot_star[:,0]-x0)**2+(snap_shot_star[:,1]-y0)**2+
               (snap_shot_star[:,2]-z0)**2)
   ids = dstar <= R0
   inlstar = snap_shot_star[ids,:]
   inlmass_star = snap_mass_star[ids]

dDM = np.sqrt((snap_shot_DM[:,0]-x0)**2+(snap_shot_DM[:,1]-y0)**2+
               (snap_shot_DM[:,2]-z0)**2)
iD = dDM <= R0
inlDM = snap_shot_DM[iD,:]
ddisk = np.sqrt((snap_shot_disk[:,0]-x0)**2+(snap_shot_disk[:,1]-y0)**2+
               (snap_shot_disk[:,2]-z0)**2)
idisk = ddisk <= R0
inldisk = snap_shot_disk[idisk,:]
dbulge = np.sqrt((snap_shot_bulge[:,0]-x0)**2+(snap_shot_bulge[:,1]-y0)**2+
               (snap_shot_bulge[:,2]-z0)**2)
ibu = dbulge <= R0
inlbulge = snap_shot_bulge[ibu,:]

if snap_N[-1] !=0:
   dbndry = np.sqrt((snap_shot_bndry[:,0]-x0)**2+(snap_shot_bndry[:,1]-y0)**2+
               (snap_shot_bndry[:,2]-z0)**2)
   idb = dbndry <= R0
   inlbndry = snap_shot_bndry[idb,:]
   inlmass_BH = snap_mass_BH[idb]
#try to find the most densitive area in the central area(which will set as BCG)
#R_range = R0/18.0
R_range = 200.0
d_gas = np.sqrt((snap_shot_gas[:,0]-x0)**2+(snap_shot_gas[:,1]-y0)**2+
               (snap_shot_gas[:,2]-z0)**2)
ig_ = d_gas <= R_range
inl_gas = snap_shot_gas[ig_,:]
if snap_N[-2] !=0:
   d_star = np.sqrt((snap_shot_star[:,0]-x0)**2+(snap_shot_star[:,1]-y0)**2+
               (snap_shot_star[:,2]-z0)**2)
   ids_ = d_star <= R_range
   inl_star = snap_shot_star[ids_,:]
##central postion and density of star
num_bins = np.ceil(R_range*2/resolution)
hist_star,edge_star = np.histogramdd(inl_star, bins = (num_bins,num_bins,num_bins))
bin_x_star = np.array(edge_star[0])
bin_y_star = np.array(edge_star[1])
bin_z_star = np.array(edge_star[2])
#to find the first five density center and then choose the closed one to the halo center
maxN = 5
cen_po_star = np.zeros((maxN,3),dtype = np.float)
hist_use = hist_star
for k in range(maxN):
    is_max = np.unravel_index(np.argmax(hist_use, axis=None), hist_use.shape)
    cen_star = hist_use[is_max]
    cenxstar = (bin_x_star[is_max[0]+1]+bin_x_star[is_max[0]])/2.0
    cenystar = (bin_y_star[is_max[1]+1]+bin_y_star[is_max[1]])/2.0
    cenzstar = (bin_z_star[is_max[2]+1]+bin_z_star[is_max[2]])/2.0
    cen_po_star[k,:] = np.array([cenxstar,cenystar,cenzstar])
    hist_use[is_max] = 0.0
compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+(cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
ismin = find.find1d(compare_d,np.min(compare_d))
cen_x_star = cen_po_star[ismin,0]
cen_y_star = cen_po_star[ismin,1]
cen_z_star = cen_po_star[ismin,2]  
# next,give radius and calculate the properties,in units kpc/h
"""
according to Milky Way,assuming the BCG radiu is about 30kpc/h,but for finding the cut point
give a 100kpc/h,to find the mean density profile
"""
R_BCG = size_BCG
r_bcg = np.logspace(-2,np.log10(R_BCG),N_size)
n_r = len(r_bcg)
inl_m_tot = np.zeros(n_r,dtype = np.float)
inl_m_g = np.zeros(n_r,dtype = np.float)
#above three mass instand:mass of a gas particle,star particle and DM particle aspectlly
r_gas = np.sqrt((inlgas[:,0]-cen_x_star)**2+(inlgas[:,1]-cen_y_star)**2+
                (inlgas[:,2]-cen_z_star)**2)
for p in range(n_r):
    ia = r_gas <= r_bcg[p]
    rgas = r_gas[ia]
    inl_m_g[p] = np.sum(inlmass_gas[ia]*10**10)
if snap_N[-2] !=0:
    inl_m_s = np.zeros(n_r,dtype = np.float)
    r_star = np.sqrt((inlstar[:,0]-cen_x_star)**2+(inlstar[:,1]-cen_y_star)**2+
                 (inlstar[:,2]-cen_z_star)**2)
    for p in range(n_r):
        ib = r_star <= r_bcg[p]
        rstar = r_star[ib]
        inl_m_s[p] = np.sum(inlmass_star[ib]*10**10)
if snap_N[-1] !=0:
    inl_BH = {}
    inl_BH_m = np.zeros(n_r,dtype = np.float)
    r_BH = np.sqrt((inlbndry[:,0]-cen_x_star)**2+(inlbndry[:,1]-cen_y_star)**2+
                 (inlbndry[:,2]-cen_z_star)**2)
    for p in range(n_r):
        ic = r_BH <=r_bcg[p]
        rbh = r_BH[ic]
        inl_BH[p] = inlbndry[ic,:]
        inl_BH_m[p] = np.sum(inlmass_BH[ic]*10**10)
if snap_N[-1] !=0 and snap_N[-2] !=0:
    inl_m_tot = inl_m_g+ inl_m_s + inl_BH_m
if snap_N[-1] ==0 and snap_N[-2] !=0:
    inl_m_tot = inl_m_g+ inl_m_s
if snap_N[-1]==0 and snap_N[-2] ==0:
    inl_m_tot = inl_m_g
inl_m_b = inl_m_g + inl_m_s
inl_rho_b = np.zeros(n_r-1,dtype = np.float)
inl_rho_g = np.zeros(n_r-1,dtype = np.float)
inl_rho_s = np.zeros(n_r-1,dtype = np.float)
for p in range(n_r-1):
    dr = r_bcg[p+1] - r_bcg[p]
    inl_rho_g[p] = (inl_m_g[p+1]-inl_m_g[p])/(4*np.pi*dr*r_bcg[p]**2)
    if snap_N[-2] !=0:
        inl_rho_s[p] = (inl_m_s[p+1]-inl_m_s[p])/(4*np.pi*dr*r_bcg[p]**2)
inl_rho_b = inl_rho_g + inl_rho_s
mean_rho_g = inl_m_g/(4*np.pi*r_bcg**3/3)
if snap_N[-2] !=0:
    mean_rho_s = inl_m_s/(4*np.pi*r_bcg**3/3)
mean_rho_b = inl_m_b/(4*np.pi*r_bcg**3/3)
def fig_center(t):
    ## figure the results of MUSIC
    if snap_N[-2] !=0:
        plt.hist2d(inlstar[:,0],inlstar[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        plt.scatter(cen_x_star,cen_y_star,s = 25, c = 'm',marker = 'o',)
        plt.scatter(x0,y0,s = 25,c = 'k',marker = '+')
        plt.xlim(x0-R0,x0+R0)
        plt.ylim(y0-R0,y0+R0)
        plt.xticks([x0-R0,x0,x0+R0])
        plt.yticks([y0,y0+R0],rotation = 90)
        plt.title('total view')
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$r-kpc/h$')
        plt.axes().set_aspect('equal')
        plt.savefig('BCG_star_mu_No1_z_128',dpi=600)
        plt.show()
    #plt.hist2d(inlDM[:,0],inlDM[:,1],bins = [1000,1000],
    #           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
    #plt.hist2d(inlgas[:,0],inlgas[:,1],bins = [1000,1000],
    #           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        #plt.hist2d(inlstar[:,0],inlstar[:,1],bins = [1000,1000],
        #       cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
        plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [num_bins,num_bins],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
    #hsc.circles(cen_po_gas[0],cen_po_gas[1],s = resolution, c = 'k')
    #hsc.circles(cen_po_star[0],cen_po_star[1],s = resolution, c = 'm')
    #hsc.circles(x0,y0,s = 5*Nr*resolution ,c = 'b',alpha = 0.35)
    A = [[cen_x_star-10*Nr*resolution,cen_x_star+10*Nr*resolution],[cen_x_star,cen_x_star]]
    A = np.array(A)
    B = [[cen_y_star,cen_y_star],[cen_y_star-10*Nr*resolution,cen_y_star+10*Nr*resolution]]
    B = np.array(B)
    plt.plot(A[0],B[0],'k-',lw = 1,alpha = 0.35,label = r'$20kpc/h_{per line}$')
    plt.plot(A[1],B[1],'k-',lw = 1,alpha = 0.35)
    C = [[x0-10*Nr*resolution,x0+10*Nr*resolution],[x0,x0]]
    D = [[y0,y0],[y0-10*Nr*resolution,y0+10*Nr*resolution]]
    plt.plot(C[0],D[0],'b-',lw = 2,alpha = 0.35,label = r'$20kpc/h_{per line}$')
    plt.plot(C[1],D[1],'b-',lw = 2,alpha = 0.35)
    if snap_N[-1] !=0:
        plt.scatter(inl_BH[n_r-1][:,0],inl_BH[-1][:,1],s = 25,c = 'k')
    plt.xlim(x0-R_range,x0+R_range)
    plt.ylim(y0-R_range,y0+R_range)
    plt.legend(loc = 4,fontsize=5)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$r-kpc/h$')
    plt.xticks([x0-R_range,x0,x0+R_range],size = 10)
    plt.yticks([y0,y0+R_range],rotation = 90,size = 10)
    plt.axes().set_aspect('equal')
    plt.title(r'$MUSIC_{%.3f}-resolution_{%.3f}$'%(snap_z,resolution))
    plt.savefig('BCG_Multi_center_mu No1_z_128',dpi=600)
    plt.show()
    return
#fig_center(t=True)
def fig_profile(t):
    plt.figure()
    plt.plot(r_bcg,inl_m_g,'g-',label = r'$M_g$')
    plt.plot(r_bcg,inl_m_s,'r-',label = r'$M_\ast$')
    plt.plot(r_bcg,inl_m_b,'b-',label = r'$M_b$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-1,np.max(r_bcg))
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$M-M_{\odot}/h$')
    plt.legend(loc = 2)
    plt.title(r'$BCG m-r MU_{%.3f}$'%resolution)
    plt.savefig('BCG mass profile MU',dpi = 600)
    plt.show()
    plt.close()
    plt.figure()
    plt.plot(r_bcg[:-1],inl_rho_g,'g-',label = r'$\rho_g$')
    plt.plot(r_bcg[:-1],inl_rho_s,'r-',label = r'$\rho_{\ast}$')
    plt.plot(r_bcg[:-1],inl_rho_b,'b-',label = r'$\rho_b$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-1,np.max(r_bcg))
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')
    plt.legend(loc = 1)
    plt.title(r'$BCG \rho-r MU_{%.3f}$'%resolution)
    plt.savefig('BCG density profile MU',dpi = 600)
    plt.show()
    plt.close()
    plt.figure()
    plt.plot(r_bcg,mean_rho_g,'g-',label = r'$\bar{\rho_g} MU$')
    plt.plot(r_bcg,mean_rho_s,'r-',label = r'$\bar{\rho_{\ast}} MU$')
    plt.plot(r_bcg,mean_rho_b,'b-',label = r'$\bar{\rho_b} MU$')
    plt.legend(loc = 1)
    plt.xlim(1e-1,np.max(r_bcg))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')
    plt.title(r'$BCG \bar{\rho}-r MU_{%.3f}$'%resolution)
    plt.savefig('BCG mean density MU',dpi = 600)
    plt.show()
    plt.close()
    return
#fig_profile(t = True)
###################################section 3: read the data of GX simution
snap_z_gx = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'redshift')
snap_z_gx = np.abs(snap_z_gx)
print('Now redshift is %.3f'%snap_z_gx)
snap_N_gx = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'npartTotal')
print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
      %(snap_N_gx[0],snap_N_gx[1],snap_N_gx[2],snap_N_gx[3],snap_N_gx[4],snap_N_gx[5]))
snap_name_gx = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'header')
#print(snap_name_gx)
#for type
snap_shot_gas_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','gas')
snap_shot_DM_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','dm')
try:
    snap_shot_star_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','star')
    snap_mass_star_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'mass','star')
except SystemExit:
    print('no star particles now')
#for respective position
snap_shot_bulge_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','bulge')
snap_shot_disk_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','disk')
try:
    snap_shot_bndry_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','bndry')
    snap_mass_BH_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'mass','bndry')
except SystemExit:
    print('no boundary particles now')
snap_mass_gas_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'mass','gas')
#for density profile of gas(there only density about gas in the simulation data)
#snap_dens_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'rho','gas') 
#the total mass distribution of snap_shot_128
main_halo_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
Rvir_gx = np.array(main_halo_gx['col12'])
xhalo_gx = np.array(main_halo_gx['col6'])
yhalo_gx = np.array(main_halo_gx['col7'])
zhalo_gx = np.array(main_halo_gx['col8'])
x0_gx = xhalo_gx[_id_]
y0_gx = yhalo_gx[_id_]
z0_gx = zhalo_gx[_id_]
R0_gx = Rvir_gx[_id_]
dgas_gx = np.sqrt((snap_shot_gas_gx[:,0]-x0_gx)**2+(snap_shot_gas_gx[:,1]-y0_gx)**2+
               (snap_shot_gas_gx[:,2]-z0_gx)**2)
ig_gx = dgas_gx <= R0_gx
inlgas_gx = snap_shot_gas_gx[ig_gx,:]
inlmass_gas_gx = snap_mass_gas_gx[ig_gx]

if snap_N_gx[-2] !=0:
   dstar_gx = np.sqrt((snap_shot_star_gx[:,0]-x0_gx)**2+(snap_shot_star_gx[:,1]-y0_gx)**2+
               (snap_shot_star_gx[:,2]-z0_gx)**2)
   ids_gx = dstar_gx <= R0_gx
   inlstar_gx = snap_shot_star_gx[ids_gx,:]
   inlmass_star_gx = snap_mass_star_gx[ids_gx]

dDM_gx = np.sqrt((snap_shot_DM_gx[:,0]-x0_gx)**2+(snap_shot_DM_gx[:,1]-y0_gx)**2+
               (snap_shot_DM_gx[:,2]-z0_gx)**2)
iD_gx = dDM_gx <= R0_gx
inlDM_gx = snap_shot_DM_gx[iD_gx,:]
ddisk_gx = np.sqrt((snap_shot_disk_gx[:,0]-x0_gx)**2+(snap_shot_disk_gx[:,1]-y0_gx)**2+
               (snap_shot_disk_gx[:,2]-z0_gx)**2)
idisk_gx = ddisk_gx <= R0_gx
inldisk_gx = snap_shot_disk_gx[idisk_gx,:]
dbulge_gx = np.sqrt((snap_shot_bulge_gx[:,0]-x0_gx)**2+(snap_shot_bulge_gx[:,1]-y0_gx)**2+
               (snap_shot_bulge_gx[:,2]-z0_gx)**2)
ibu_gx = dbulge_gx <= R0_gx
inlbulge_gx = snap_shot_bulge_gx[ibu_gx,:]
   
if snap_N_gx[-1] !=0:
   dbndry_gx = np.sqrt((snap_shot_bndry_gx[:,0]-x0_gx)**2+(snap_shot_bndry_gx[:,1]-y0_gx)**2+
               (snap_shot_bndry_gx[:,2]-z0_gx)**2)
   idb_gx = dbndry_gx <= R0_gx
   inlbndry_gx = snap_shot_bndry_gx[idb_gx,:]
   inlmass_BH_gx = snap_mass_BH_gx[idb_gx]
#try to find the most densitive area in the central area(which will set as BCG)
#R_range_gx = R0_gx/18.0
R_range_gx = 200.0
d_gas_gx = np.sqrt((snap_shot_gas_gx[:,0]-x0_gx)**2+(snap_shot_gas_gx[:,1]-y0_gx)**2+
               (snap_shot_gas_gx[:,2]-z0_gx)**2)
i_g_gx = d_gas_gx <= R_range_gx
inl_gas_gx = snap_shot_gas_gx[i_g_gx,:]
if snap_N_gx[-2] !=0:
   d_star_gx = np.sqrt((snap_shot_star_gx[:,0]-x0_gx)**2+(snap_shot_star_gx[:,1]-y0_gx)**2+
               (snap_shot_star_gx[:,2]-z0_gx)**2)
   id_s_gx = d_star_gx <= R_range_gx
   inl_star_gx = snap_shot_star_gx[id_s_gx,:]
##central postion and density of star
num_bins_gx = np.ceil(R_range_gx*2/resolution)
hist_star_gx,edge_star_gx = np.histogramdd(inl_star_gx, bins = (num_bins_gx,num_bins_gx,num_bins_gx))
bin_x_star_gx = np.array(edge_star_gx[0])
bin_y_star_gx = np.array(edge_star_gx[1])
bin_z_star_gx = np.array(edge_star_gx[2])
#next to find the center of BCG
maxN_gx = 5
cen_po_star_gx = np.zeros((maxN_gx,3),dtype = np.float)
hist_use_gx = hist_star_gx
for k in range(maxN_gx):
    is_max_gx = np.unravel_index(np.argmax(hist_use_gx, axis=None), hist_use_gx.shape)
    cen_star_gx = hist_use_gx[is_max_gx]
    cenxstar_gx = (bin_x_star_gx[is_max_gx[0]+1]+bin_x_star_gx[is_max_gx[0]])/2.0
    cenystar_gx = (bin_y_star_gx[is_max_gx[1]+1]+bin_y_star_gx[is_max_gx[1]])/2.0
    cenzstar_gx = (bin_z_star_gx[is_max_gx[2]+1]+bin_z_star_gx[is_max_gx[2]])/2.0
    cen_po_star_gx[k,:] = np.array([cenxstar_gx,cenystar_gx,cenzstar_gx])
    hist_use_gx[is_max_gx] = 0.0
compare_d_gx= np.sqrt((cen_po_star_gx[:,0]-x0_gx)**2+
                      (cen_po_star_gx[:,1]-y0_gx)**2+(cen_po_star_gx[:,2]-z0_gx)**2)
ismin_gx = find.find1d(compare_d_gx,np.min(compare_d_gx))
cen_x_star_gx = cen_po_star_gx[ismin_gx,0]
cen_y_star_gx = cen_po_star_gx[ismin_gx,1]
cen_z_star_gx = cen_po_star_gx[ismin_gx,2]  
'''
##central postion and density of gas
hist_gas_gx,edge_gas_gx = np.histogramdd(inl_gas_gx, bins = (num_bins_gx,num_bins_gx,num_bins_gx))
bin_x_gas_gx = np.array(edge_gas_gx[0])
bin_y_gas_gx = np.array(edge_gas_gx[1])
bin_z_gas_gx = np.array(edge_gas_gx[2])
ig_max_gx = np.unravel_index(np.argmax(hist_gas_gx, axis=None), hist_gas_gx.shape)
cen_gas_gx = hist_gas_gx[ig_max_gx]
cen_x_gas_gx = (bin_x_gas_gx[ig_max_gx[0]+1]+bin_x_gas_gx[ig_max_gx[0]])/2.0
cen_y_gas_gx = (bin_y_gas_gx[ig_max_gx[1]+1]+bin_y_gas_gx[ig_max_gx[1]])/2.0
cen_z_gas_gx = (bin_z_gas_gx[ig_max_gx[2]+1]+bin_z_gas_gx[ig_max_gx[2]])/2.0
cen_po_gas_gx = np.array([cen_x_gas_gx,cen_y_gas_gx,cen_z_gas_gx])
'''
# next,give radius and calculate the properties,in units kpc/h
"""
according to Milky Way,assuming the BCG radiu is about 30kpc/h,but for finding the cut point
give a 100kpc/h,to find the mean density profile
"""
R_BCG_gx = size_BCG
r_bcg_gx = np.logspace(-2,np.log10(R_BCG_gx),N_size)
n_r_gx = len(r_bcg_gx)
inl_m_tot_gx = np.zeros(n_r_gx,dtype = np.float)
inl_m_g_gx = np.zeros(n_r_gx,dtype = np.float)
#above three mass instand:mass of a gas particle,star particle and DM particle aspectlly
r_gas_gx = np.sqrt((inlgas_gx[:,0]-cen_x_star_gx)**2+(inlgas_gx[:,1]-cen_y_star_gx)**2+
                (inlgas_gx[:,2]-cen_z_star_gx)**2)
for q in range(n_r_gx):
    ia_gx = r_gas_gx <= r_bcg_gx[q]
    rgas_gx = r_gas_gx[ia_gx]
    inl_m_g_gx[q] = np.sum(inlmass_gas_gx[ia_gx]*10**10)
if snap_N_gx[-2] !=0:
    inl_m_s_gx = np.zeros(n_r_gx,dtype = np.float)
    r_star_gx = np.sqrt((inlstar_gx[:,0]-cen_x_star_gx)**2+(inlstar_gx[:,1]-cen_y_star_gx)**2+
                 (inlstar_gx[:,2]-cen_z_star_gx)**2)
    for q in range(n_r_gx):
        ib_gx = r_star_gx <= r_bcg_gx[q]
        rstar_gx = r_star_gx[ib_gx]
        inl_m_s_gx[q] = np.sum(inlmass_star_gx[ib_gx]*10**10)
if snap_N_gx[-1] !=0:
    inl_BH_gx = {}
    inl_BH_m_gx = np.zeros(n_r_gx,dtype = np.float)
    r_BH_gx = np.sqrt((inlbndry_gx[:,0]-cen_x_star_gx)**2+(inlbndry_gx[:,1]-cen_y_star_gx)**2+
                 (inlbndry_gx[:,2]-cen_z_star_gx)**2)
    for q in range(n_r_gx):
        ic_gx = r_BH_gx <=r_bcg_gx[q]
        rbh_gx = r_BH_gx[ic_gx]
        inl_BH_gx[q] = inlbndry_gx[ic_gx,:]
        inl_BH_m_gx[q] = np.sum(inlmass_BH_gx[ic_gx]*10**10)
if snap_N_gx[-1] !=0 and snap_N_gx[-2] !=0:
    inl_m_tot_gx = inl_m_g_gx + inl_m_s_gx + inl_BH_m_gx
if snap_N_gx[-1] ==0 and snap_N_gx[-2] !=0:
    inl_m_tot_gx = inl_m_g_gx + inl_m_s_gx
if snap_N_gx[-1]==0 and snap_N_gx[-2] ==0:
    inl_m_tot_gx = inl_m_g_gx
inl_m_b_gx = inl_m_g_gx + inl_m_s_gx
inl_rho_b_gx = np.zeros(n_r_gx-1,dtype = np.float)
inl_rho_g_gx = np.zeros(n_r_gx-1,dtype = np.float)
inl_rho_s_gx = np.zeros(n_r_gx-1,dtype = np.float)
for q in range(n_r_gx-1):
    dr_gx = r_bcg_gx[q+1] - r_bcg_gx[q]
    inl_rho_g_gx[q] = (inl_m_g_gx[q+1]-inl_m_g_gx[q])/(4*np.pi*dr_gx*r_bcg_gx[q]**2)
    if snap_N_gx[-2] !=0:
        inl_rho_s_gx[q] = (inl_m_s_gx[q+1]-inl_m_s_gx[q])/(4*np.pi*dr_gx*r_bcg_gx[q]**2)
inl_rho_b_gx = inl_rho_g_gx + inl_rho_s_gx
mean_rho_g_gx = inl_m_g_gx/(4*np.pi*r_bcg_gx**3/3)
if snap_N_gx[-2] !=0:
    mean_rho_s_gx = inl_m_s_gx/(4*np.pi*r_bcg_gx**3/3)
mean_rho_b_gx = inl_m_b_gx/(4*np.pi*r_bcg_gx**3/3)
def fig_center_gx(p):
    ## figure the results of GX
    if snap_N_gx[-2] !=0:
        plt.hist2d(inlstar_gx[:,0],inlstar_gx[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        plt.scatter(cen_x_star_gx,cen_y_star_gx,s = 25, c = 'm',marker = 'o',)
        plt.scatter(x0_gx,y0_gx,s = 25,c = 'k',marker = '+')
        plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
        plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
        plt.yticks([y0_gx-R0_gx,y0_gx,y0_gx+R0_gx],rotation = 90)
        plt.xlabel(r'$r-kpc/h$')
        plt.ylabel(r'$r-kpc/h$')
        plt.axes().set_aspect('equal')
        plt.savefig('BCG_star_gx_No1_z_128',dpi=600)
        plt.show()
    #plt.hist2d(inlDM_gx[:,0],inlDM_gx[:,1],bins = [1000,1000],
    #           cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
    #plt.hist2d(inlgas_gx[:,0],inlgas_gx[:,1],bins = [1000,1000],
    #           cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
    if snap_N_gx[-2] !=0:
        #plt.hist2d(inlstar_gx[:,0],inlstar_gx[:,1],bins = [1000,1000],
        #       cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
        plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [num_bins_gx,num_bins_gx],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
    #hsc.circles(cen_po_gas_gx[0],cen_po_gas_gx[1],s = resolution, c = 'k')
    #hsc.circles(cen_po_star_gx[0],cen_po_star_gx[1],s = resolution, c = 'm')
    #hsc.circles(x0_gx,y0_gx,s = 5*Nr*resolution,c = 'b',alpha = 0.35)
    A_gx = [[cen_x_star_gx-10*Nr*resolution,cen_x_star_gx+10*Nr*resolution],[cen_x_star_gx,cen_x_star_gx]]
    A_gx = np.array(A_gx)
    B_gx = [[cen_y_star_gx,cen_y_star_gx],[cen_y_star_gx-10*Nr*resolution,cen_y_star_gx+10*Nr*resolution]]
    B_gx = np.array(B_gx)
    plt.plot(A_gx[0],B_gx[0],'k-',lw = 1,alpha = 0.35,label = r'$20kpc/h_{per line}$')
    plt.plot(A_gx[1],B_gx[1],'k-',lw = 1,alpha = 0.35)
    C_gx = [[x0_gx-10*Nr*resolution,x0_gx+10*Nr*resolution],[x0_gx,x0_gx]]
    D_gx = [[y0_gx,y0_gx],[y0_gx-10*Nr*resolution,y0_gx+10*Nr*resolution]]
    plt.plot(C_gx[0],D_gx[0],'b-',lw = 2,alpha = 0.35,label = r'$20kpc/h_{per line}$')
    plt.plot(C_gx[1],D_gx[1],'b-',lw = 2,alpha = 0.35)
    if snap_N_gx[-1] !=0:
        plt.scatter(inl_BH_gx[n_r_gx-1][:,0],inl_BH_gx[n_r_gx-1][:,1],s = 10,c = 'k')
    plt.xlim(x0_gx-R_range_gx ,x0_gx+R_range_gx )
    plt.ylim(y0_gx-R_range_gx ,y0_gx+R_range_gx )
    plt.legend(loc = 4,fontsize=5)
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$r-kpc/h$')
    plt.xticks([x0_gx-R_range_gx,x0_gx,x0_gx+R_range_gx],size = 10)
    plt.yticks([y0_gx,y0_gx+R_range_gx],rotation = 90,size = 10)
    plt.axes().set_aspect('equal')
    plt.title(r'$GadgetX_{%.3f}-resolution_{%.3f}$'%(snap_z_gx,resolution))
    plt.savefig('BCG_Multi_center_gx_No1_z_128',dpi=600)
    plt.show()
    return
#fig_center_gx(p = True)
def fig_profile_gx(p):
    plt.plot(r_bcg_gx,inl_m_g_gx,'g-',label = r'$M_g$')
    plt.plot(r_bcg_gx,inl_m_s_gx,'r-',label = r'$M_\ast$')
    plt.plot(r_bcg_gx,inl_m_b_gx,'b-',label = r'$M_b$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-1,np.max(r_bcg_gx))
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$M-M_{\odot}/h$')
    plt.legend(loc = 2)
    plt.title(r'$BCG m-r GX_{%.3f}$'%resolution)
    plt.savefig('BCG mass profile GX',dpi = 600)
    plt.show()
    plt.plot(r_bcg_gx[:-1],inl_rho_g_gx,'g-',label = r'$\rho_g$')
    plt.plot(r_bcg_gx[:-1],inl_rho_s_gx,'r-',label = r'$\rho_{\ast}$')
    plt.plot(r_bcg_gx[:-1],inl_rho_b_gx,'b-',label = r'$\rho_b$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-1,np.max(r_bcg_gx))
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')
    plt.legend(loc = 1)
    plt.title(r'$BCG \rho-r GX_{%.3f}$'%resolution)
    plt.savefig('BCG density profile GX',dpi = 600)
    plt.show()
    plt.plot(r_bcg,mean_rho_g_gx,'g-',label = r'$\bar{\rho_g}$')
    plt.plot(r_bcg,mean_rho_s_gx,'r-',label = r'$\bar{\rho_{\ast}}$')
    plt.plot(r_bcg,mean_rho_b_gx,'b-',label = r'$\bar{\rho_b}$')
    plt.legend(loc = 1)
    plt.xlim(1e-1,np.max(r_bcg))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')
    plt.title(r'$BCG \bar{\rho}-r GX_{%.3f}$'%resolution)
    plt.savefig('BCG mean density GX',dpi = 600)
    plt.show()
    return
#fig_profile_gx(p = True)
def comparation(q):
    dd = 128
    plt.figure()
    plt.plot(r_bcg,inl_m_g,'g-',label = r'$M_g MU$')
    plt.plot(r_bcg,inl_m_s,'r-',label = r'$M_\ast MU$')
    plt.plot(r_bcg,inl_m_b,'b-',label = r'$M_b MU$')
    plt.plot(r_bcg_gx,inl_m_g_gx,'g--',label = r'$M_g GX$')
    plt.plot(r_bcg_gx,inl_m_s_gx,'r--',label = r'$M_\ast GX$')
    plt.plot(r_bcg_gx,inl_m_b_gx,'b--',label = r'$M_b GX$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-1,np.max(r_bcg))
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$M-M_{\odot}/h$')
    plt.legend(loc = 2)
    plt.title(r'$BCG m-r Re_{%.3f}$'%resolution)
    plt.savefig('BCG m-comparartion z_%.0f.png'%dd, dpi = 600)
    plt.show()
    plt.close()
    plt.figure()
    plt.plot(r_bcg[:-1],inl_rho_g,'g-',label = r'$\rho_g MU$')
    plt.plot(r_bcg[:-1],inl_rho_s,'r-',label = r'$\rho_{\ast} MU$')
    plt.plot(r_bcg[:-1],inl_rho_b,'b-',label = r'$\rho_b MU$')
    plt.plot(r_bcg_gx[:-1],inl_rho_g_gx,'g--',label = r'$\rho_g GX$')
    plt.plot(r_bcg_gx[:-1],inl_rho_s_gx,'r--',label = r'$\rho_{\ast} GX$')
    plt.plot(r_bcg_gx[:-1],inl_rho_b_gx,'b--',label = r'$\rho_b GX$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-1,np.max(r_bcg))
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$\rho-[M_{\odot} h^2/{kpc^3}]$')
    plt.legend(loc = 1)
    plt.title(r'$BCG \rho-r Re_{%.3f}$'%resolution)
    plt.savefig('BCG rho-comparation z_%.0f.png'%dd, dpi = 600)
    plt.show()
    plt.close()
    plt.figure()
    plt.plot(r_bcg,mean_rho_g,'g-',label = r'$\bar{\rho_g} MU$')
    plt.plot(r_bcg,mean_rho_s,'r-',label = r'$\bar{\rho_{\ast}} MU$')
    plt.plot(r_bcg,mean_rho_b,'b-',label = r'$\bar{\rho_b} MU$')
    plt.plot(r_bcg,mean_rho_g_gx,'g--',label = r'$\bar{\rho_g} GX$')
    plt.plot(r_bcg,mean_rho_s_gx,'r--',label = r'$\bar{\rho_{\ast}} GX$')
    plt.plot(r_bcg,mean_rho_b_gx,'b--',label = r'$\bar{\rho_b} GX$')
    plt.legend(loc = 1)
    plt.xlim(1e-1,np.max(r_bcg))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$r-kpc/h$')
    plt.ylabel(r'$\bar{\rho} [M_{\odot} h^2/kpc^3]$')
    plt.title(r'$BCG \bar{\rho}-r Re_{%.3f}$'%resolution)
    plt.savefig('BCG mean density comparation',dpi = 600)
    plt.show()
    plt.close()
    return
#comparation(q = True)