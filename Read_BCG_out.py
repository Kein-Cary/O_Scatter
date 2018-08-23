#this file try to find the BCG in a given halo(to find the most dengsity part of given halo)
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
##two files save the snap_shot data
#'D:/mask/snapshot/MUSIC/snap_000'
#'D:/mask/snapshot/GX/snap_000'
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
No_snap = '128'
#######################section1:read out the position data of MUSIC
snap_z = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'redshift')
snap_z = np.abs(snap_z)
print('Now redshift is %.3f'%snap_z)
snap_N = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'npartTotal')
print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
      %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
snap_name = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'header')
#for type
snap_shot_gas = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','gas')
snap_shot_DM = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','dm')
try:
    snap_shot_star = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','star')
except SystemExit:
    print('no star particles now')
#for respective position
snap_shot_bulge = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','bulge')
snap_shot_disk = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','disk')
try:
    snap_shot_bndry = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'pos','bndry')
except SystemExit:
    print('no boundary particles now')
#for density profile of gas(there only density about gas in the simulation data)
#snap_dens = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'rho','gas') 
#the total mass distribution of snap_shot_128
##next,try to find the particles belong to halo 128000000000001(in MUSIC simulation)
main_halo = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
Rvir = np.array(main_halo['col12'])
xhalo = np.array(main_halo['col6'])
yhalo = np.array(main_halo['col7'])
zhalo = np.array(main_halo['col8'])
x0 = xhalo[0]
y0 = yhalo[0]
z0 = zhalo[0]
R0 = Rvir[0]
dgas = np.sqrt((snap_shot_gas[:,0]-x0)**2+(snap_shot_gas[:,1]-y0)**2+
               (snap_shot_gas[:,2]-z0)**2)
ig = dgas <= R0
inl_gas = snap_shot_gas[ig,:]

dDM = np.sqrt((snap_shot_DM[:,0]-x0)**2+(snap_shot_DM[:,1]-y0)**2+
               (snap_shot_DM[:,2]-z0)**2)
iD = dDM <= R0
inl_DM = snap_shot_DM[iD,:]

ddisk = np.sqrt((snap_shot_disk[:,0]-x0)**2+(snap_shot_disk[:,1]-y0)**2+
               (snap_shot_disk[:,2]-z0)**2)
idisk = ddisk <= R0
inl_disk = snap_shot_disk[idisk,:]

dbulge = np.sqrt((snap_shot_bulge[:,0]-x0)**2+(snap_shot_bulge[:,1]-y0)**2+
               (snap_shot_bulge[:,2]-z0)**2)
ibu = dbulge <= R0
inl_bulge = snap_shot_bulge[ibu,:]

if snap_N[-2] !=0:
   dstar = np.sqrt((snap_shot_star[:,0]-x0)**2+(snap_shot_star[:,1]-y0)**2+
               (snap_shot_star[:,2]-z0)**2)
   ids = dstar <= R0
   inl_star = snap_shot_star[ids,:]
#try to find the most massive aera,inclode number density and coordinate
##central postion and density of gas
hist_gas,edge_gas = np.histogramdd(inl_gas, bins = (50,50,50))
bin_x_gas = np.array(edge_gas[0])
bin_y_gas = np.array(edge_gas[1])
bin_z_gas = np.array(edge_gas[2])
cen_gas = np.zeros((inl_gas.shape[0],3),dtype = np.float)
N_gas = len(hist_gas[hist_gas!=0])
dens_gas = np.zeros(N_gas,dtype = np.float)
hist_1 = hist_gas
for k in range(N_gas):
    idd = np.unravel_index(np.argmax(hist_1, axis=None), hist_1.shape)
    dens_gas[k] = hist_gas[idd]
    idd = np.array(idd)
    if idd[0] ==100 :
        cen_gas[k,0] = bin_x_gas[idd[0]]
    else:
        cen_gas[k,0] = (bin_x_gas[idd[0]+1]+bin_x_gas[idd[0]])/2
    if idd[1] ==100 :
        cen_gas[k,1] = bin_y_gas[idd[1]]
    else:
        cen_gas[k,1] = (bin_y_gas[idd[1]+1]+bin_y_gas[idd[1]])/2   
    if idd[2] ==100 :
        cen_gas[k,2] = bin_z_gas[idd[2]]
    else:
        cen_gas[k,2] = (bin_z_gas[idd[2]+1]+bin_z_gas[idd[2]])/2  
    hist_1[idd[0],idd[1],idd[2]] = 0
cen_po_gas = cen_gas[cen_gas!=0]
cen_po_gas = cen_po_gas.reshape(np.int(len(cen_po_gas)/3),3)
with h5py.File('cen_po_gas_music.h5','w') as f:
    f['a'] = np.array(cen_po_gas)
with h5py.File('cen_po_gas_music.h5') as f:
    for k in range(np.int(len(cen_po_gas)/3)):
        f['a'][k,:] = cen_po_gas[k,:]
##central postion and density of star
hist_star,edge_star = np.histogramdd(inl_star, bins = (50,50,50))
bin_x_star = np.array(edge_star[0])
bin_y_star = np.array(edge_star[1])
bin_z_star = np.array(edge_star[2])
cen_star = np.zeros((inl_star.shape[0],3),dtype = np.float)
N_star = len(hist_star[hist_star!=0])
dens_star = np.zeros(N_star,dtype = np.float)
hist_2 = hist_star
for k in range(N_star):
    idd_star = np.unravel_index(np.argmax(hist_2, axis=None), hist_2.shape)
    dens_star[k] = hist_star[idd_star]
    idd_star = np.array(idd_star)
    if idd_star[0] ==100 :
        cen_star[k,0] = bin_x_star[idd_star[0]]
    else:
        cen_star[k,0] = (bin_x_star[idd_star[0]+1]+bin_x_star[idd_star[0]])/2
    if idd_star[1] ==100 :
        cen_star[k,1] = bin_y_star[idd_star[1]]
    else:
        cen_star[k,1] = (bin_y_star[idd_star[1]+1]+bin_y_star[idd_star[1]])/2   
    if idd_star[2] ==100 :
        cen_star[k,2] = bin_z_star[idd_star[2]]
    else:
        cen_star[k,2] = (bin_z_star[idd_star[2]+1]+bin_z_star[idd_star[2]])/2  
    hist_2[idd_star[0],idd_star[1],idd_star[2]] = 0
cen_po_star = cen_star[cen_star!=0]
cen_po_star = cen_po_star.reshape(np.int(len(cen_po_star)/3),3)
with h5py.File('cen_po_star_music.h5','w') as f:
    f['a'] = np.array(cen_po_star)
with h5py.File('cen_po_star_music.h5') as f:
    for k in range(np.int(len(cen_po_star)/3)):
        f['a'][k,:] = cen_po_star[k,:]
with h5py.File('music_gas.h5','w') as f:
    f['a'] = np.array(dens_gas)
with h5py.File('music_star.h5','w') as f:
    f['a'] = np.array(dens_star) 
'''
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(cen_po_star[:,0],cen_po_star[:,1],cen_po_star[:,2],s=0.5,c = dens_star,
           marker='o',cmap='plasma',vmin=np.min(dens_star),vmax=np.max(dens_star),
           norm = mpl.colors.LogNorm(),alpha=0.2)
plt.savefig('star_distribution',dpi=600)
'''
###################################section 3: read the data of GX simution
snap_z_gx = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'redshift')
snap_z_gx = np.abs(snap_z_gx)
print('Now redshift is %.3f'%snap_z_gx)
snap_N_gx = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'npartTotal')
print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
      %(snap_N_gx[0],snap_N_gx[1],snap_N_gx[2],snap_N_gx[3],snap_N_gx[4],snap_N_gx[5]))
snap_name_gx = pygdr.readheader('D:/mask/snapshot/GX/snap_%s'%No_snap,'header')
#for type
snap_shot_gas_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','gas')
snap_shot_DM_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','dm')
try:
    snap_shot_star_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','star')
except SystemExit:
    print('no star particles now')
#for respective position
snap_shot_bulge_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','bulge')
snap_shot_disk_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','disk')
try:
    snap_shot_bndry_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'pos','bndry')
except SystemExit:
    print('no boundary particles now')
#for density profile of gas(there only density about gas in the simulation data)
#snap_dens_gx = pygdr.readsnap('D:/mask/snapshot/GX/snap_%s'%No_snap,'rho','gas') 
#the total mass distribution of snap_shot_128
##next,try to find the particles belong to halo 128000000000001(in GX simulation)
main_halo_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
Rvir_gx = np.array(main_halo_gx['col12'])
xhalo_gx = np.array(main_halo_gx['col6'])
yhalo_gx = np.array(main_halo_gx['col7'])
zhalo_gx = np.array(main_halo_gx['col8'])
x0_gx = xhalo_gx[0]
y0_gx = yhalo_gx[0]
z0_gx = zhalo_gx[0]
R0_gx = Rvir_gx[0]
dgas_gx = np.sqrt((snap_shot_gas_gx[:,0]-x0_gx)**2+(snap_shot_gas_gx[:,1]-y0_gx)**2+
               (snap_shot_gas_gx[:,2]-z0_gx)**2)
ig_gx = dgas_gx <= R0_gx
inl_gas_gx = snap_shot_gas_gx[ig_gx,:]

dDM_gx = np.sqrt((snap_shot_DM_gx[:,0]-x0_gx)**2+(snap_shot_DM_gx[:,1]-y0_gx)**2+
               (snap_shot_DM_gx[:,2]-z0_gx)**2)
iD_gx = dDM_gx <= R0_gx
inl_DM_gx = snap_shot_DM_gx[iD_gx,:]

ddisk_gx = np.sqrt((snap_shot_disk_gx[:,0]-x0_gx)**2+(snap_shot_disk_gx[:,1]-y0_gx)**2+
               (snap_shot_disk_gx[:,2]-z0_gx)**2)
idisk_gx = ddisk_gx <= R0_gx
inl_disk_gx = snap_shot_disk_gx[idisk_gx,:]

dbulge_gx = np.sqrt((snap_shot_bulge_gx[:,0]-x0_gx)**2+(snap_shot_bulge_gx[:,1]-y0_gx)**2+
               (snap_shot_bulge_gx[:,2]-z0_gx)**2)
ibu_gx = dbulge_gx <= R0_gx
inl_bulge_gx = snap_shot_bulge_gx[ibu_gx,:]

if snap_N_gx[-2] !=0:
   dstar_gx = np.sqrt((snap_shot_star_gx[:,0]-x0_gx)**2+(snap_shot_star_gx[:,1]-y0_gx)**2+
               (snap_shot_star_gx[:,2]-z0_gx)**2)
   ids_gx = dstar_gx <= R0_gx
   inl_star_gx = snap_shot_star_gx[ids_gx,:]
#try to find the most massive aera,inclode number density and coordinate
##central postion and density of gas
hist_gas_gx,edge_gas_gx = np.histogramdd(inl_gas_gx, bins = (50,50,50))
bin_x_gas_gx = np.array(edge_gas_gx[0])
bin_y_gas_gx = np.array(edge_gas_gx[1])
bin_z_gas_gx = np.array(edge_gas_gx[2])
cen_gas_gx = np.zeros((inl_gas_gx.shape[0],3),dtype = np.float)
N_gas_gx = len(hist_gas_gx[hist_gas_gx!=0])
dens_gas_gx = np.zeros(N_gas_gx,dtype = np.float)
hist_1_gx = hist_gas_gx
for k in range(N_gas_gx):
    idd_gx = np.unravel_index(np.argmax(hist_1_gx, axis=None), hist_1_gx.shape)
    dens_gas_gx[k] = hist_gas_gx[idd_gx]
    idd_gx = np.array(idd_gx)
    if idd_gx[0] ==100 :
        cen_gas_gx[k,0] = bin_x_gas_gx[idd_gx[0]]
    else:
        cen_gas_gx[k,0] = (bin_x_gas_gx[idd_gx[0]+1]+bin_x_gas_gx[idd_gx[0]])/2
    if idd_gx[1] ==100 :
        cen_gas_gx[k,1] = bin_y_gas_gx[idd_gx[1]]
    else:
        cen_gas_gx[k,1] = (bin_y_gas_gx[idd_gx[1]+1]+bin_y_gas_gx[idd_gx[1]])/2   
    if idd_gx[2] ==100 :
        cen_gas_gx[k,2] = bin_z_gas_gx[idd_gx[2]]
    else:
        cen_gas_gx[k,2] = (bin_z_gas_gx[idd_gx[2]+1]+bin_z_gas_gx[idd_gx[2]])/2  
    hist_1_gx[idd_gx[0],idd_gx[1],idd_gx[2]] = 0
cen_po_gas_gx = cen_gas_gx[cen_gas_gx!=0]
cen_po_gas_gx = cen_po_gas_gx.reshape(np.int(len(cen_po_gas_gx)/3),3)
with h5py.File('cen_po_gas_gx.h5','w') as f:
    f['a'] = np.array(cen_po_gas_gx)
with h5py.File('cen_po_gas_gx.h5') as f:
    for k in range(np.int(len(cen_po_gas_gx)/3)):
        f['a'][k,:] = cen_po_gas_gx[k,:]
##central postion and density of star
hist_star_gx,edge_star_gx = np.histogramdd(inl_star_gx, bins = (50,50,50))
bin_x_star_gx = np.array(edge_star_gx[0])
bin_y_star_gx = np.array(edge_star_gx[1])
bin_z_star_gx = np.array(edge_star_gx[2])
cen_star_gx = np.zeros((inl_star_gx.shape[0],3),dtype = np.float)
N_star_gx = len(hist_star_gx[hist_star_gx!=0])
dens_star_gx = np.zeros(N_star_gx,dtype = np.float)
hist_2_gx = hist_star_gx
for k in range(N_star_gx):
    idd_star_gx = np.unravel_index(np.argmax(hist_2_gx, axis=None), hist_2_gx.shape)
    dens_star_gx[k] = hist_star_gx[idd_star_gx]
    idd_star_gx = np.array(idd_star_gx)
    if idd_star_gx[0] ==100 :
        cen_star_gx[k,0] = bin_x_star_gx[idd_star_gx[0]]
    else:
        cen_star_gx[k,0] = (bin_x_star_gx[idd_star_gx[0]+1]+bin_x_star_gx[idd_star_gx[0]])/2
    if idd_star_gx[1] ==100 :
        cen_star_gx[k,1] = bin_y_star_gx[idd_star_gx[1]]
    else:
        cen_star_gx[k,1] = (bin_y_star_gx[idd_star_gx[1]+1]+bin_y_star_gx[idd_star_gx[1]])/2   
    if idd_star_gx[2] ==100 :
        cen_star_gx[k,2] = bin_z_star_gx[idd_star_gx[2]]
    else:
        cen_star_gx[k,2] = (bin_z_star_gx[idd_star_gx[2]+1]+bin_z_star_gx[idd_star_gx[2]])/2  
    hist_2_gx[idd_star_gx[0],idd_star_gx[1],idd_star_gx[2]] = 0
cen_po_star_gx = cen_star_gx[cen_star_gx!=0]
cen_po_star_gx = cen_po_star_gx.reshape(np.int(len(cen_po_star_gx)/3),3)
with h5py.File('cen_po_star_gx.h5','w') as f:
    f['a'] = np.array(cen_po_star_gx)
with h5py.File('cen_po_star_gx.h5') as f:
    for k in range(np.int(len(cen_po_star_gx)/3)):
        f['a'][k,:] = cen_po_star_gx[k,:]
with h5py.File('GX_gas.h5','w') as f:
    f['a'] = np.array(dens_gas_gx)
with h5py.File('GX_star.h5','w') as f:
    f['a'] = np.array(dens_star_gx)  
'''
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax_gx = Axes3D(fig)
ax_gx.scatter(cen_po_star_gx[:,0],cen_po_star_gx[:,1],cen_po_star_gx[:,2],s=0.5,c = dens_star_gx,
           marker='o',cmap='plasma',vmin=np.min(dens_star_gx),vmax=np.max(dens_star_gx),
           norm = mpl.colors.LogNorm(),alpha=0.2)
plt.savefig('star_distribution_gx',dpi=600)
''' 