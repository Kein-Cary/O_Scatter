#this file try to point out the enter of halo,BCG,satllite galaxy
from pygadgetreader import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from handy import scatter as hsc
import numpy as np
import h5py
import pygadgetreader as pygdr
import pandas as pd
import astropy.io.ascii as asc
##two files save the snap_shot data
#'D:/mask/snapshot/MUSIC/snap_000'
#'D:/mask/snapshot/GX/snap_000'
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
##### The 1st part,read the snapshot data,and the halo information
No_snap = '128'
_id_ = 0
##set the goal snapshot and halo id
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
#for density profile of gas
#snap_dens = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'rho','gas') 
#the total mass distribution of snap_shot_128,about gas
main_halo = asc.read('D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
Rvir = np.array(main_halo['col12'])
xhalo = np.array(main_halo['col6'])
yhalo = np.array(main_halo['col7'])
zhalo = np.array(main_halo['col8'])
#try to find the subhalo belong to the first halo
halo_id = np.array(main_halo['col1'])
halo_id = np.int64(halo_id)
host_id = np.array(main_halo['col2'])
host_id = np.int64(host_id)
goal_halo = halo_id[_id_]
ix = host_id == goal_halo
po_halo = np.zeros((len(halo_id),3),dtype = np.float)
for k in range(len(ix)):
    if k == _id_ :
        po_halo[k,:] = np.array([xhalo[k],yhalo[k],zhalo[k]])
    elif ix[k] == True :
        po_halo[k,:] = np.array([xhalo[k],yhalo[k],zhalo[k]])
    else:
        po_halo[k,:] = np.array([np.inf,np.inf,np.inf]) 
po_halo = po_halo[po_halo != np.inf]
po_halo = po_halo.reshape(np.int(len(po_halo)/3),3)
#try to find the particles belong to this halo
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
if snap_N[-1] !=0:
   dbnd = np.sqrt((snap_shot_bndry[:,0]-x0)**2+(snap_shot_bndry[:,1]-y0)**2+
               (snap_shot_bndry[:,2]-z0)**2)
   ibnd = dbnd <= R0
   inl_bndry = snap_shot_bndry[ibnd,:]
##### The 2nd part,read the snapshot data,and the halo information about GX simulation
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
#for density profile of gas
#snap_dens = pygdr.readsnap('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'rho','gas') 
#the total mass distribution of snap_shot_128
##next,try to find the particles belong to halo 128000000000001(in GX simulation)
main_halo_gx = asc.read('D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos',
                  converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
Rvir_gx = np.array(main_halo_gx['col12'])
xhalo_gx = np.array(main_halo_gx['col6'])
yhalo_gx = np.array(main_halo_gx['col7'])
zhalo_gx = np.array(main_halo_gx['col8'])
#try to find the subhalo belong to the first halo
halo_id_gx = np.array(main_halo_gx['col1'])
halo_id_gx = np.int64(halo_id_gx)
host_id_gx = np.array(main_halo_gx['col2'])
host_id_gx = np.int64(host_id_gx)
goal_halo_gx = halo_id_gx[_id_]
iy = host_id_gx == goal_halo_gx
po_halo_gx = np.zeros((len(halo_id_gx),3),dtype = np.float)
for k in range(len(iy)):
    if k == _id_ :
        po_halo_gx[k,:] = np.array([xhalo_gx[k],yhalo_gx[k],zhalo_gx[k]])
    elif iy[k] == True :
        po_halo_gx[k,:] = np.array([xhalo_gx[k],yhalo_gx[k],zhalo_gx[k]])
    else:
        po_halo_gx[k,:] = np.array([np.inf,np.inf,np.inf]) 
po_halo_gx = po_halo_gx[po_halo_gx != np.inf]
po_halo_gx = po_halo_gx.reshape(np.int(len(po_halo_gx)/3),3)
#try to find the particles belong to this halo
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
if snap_N_gx[-1] !=0:
   dbnd_gx = np.sqrt((snap_shot_bndry_gx[:,0]-x0_gx)**2+(snap_shot_bndry_gx[:,1]-y0_gx)**2+
               (snap_shot_bndry_gx[:,2]-z0_gx)**2)
   ibnd_gx = dbnd_gx <= R0_gx
   inl_bndry_gx = snap_shot_bndry_gx[ibnd_gx,:]
##### The 3rd part,read the center of the particles distribution
with h5py.File('cen_po_star_music.h5') as f:
    cen_po_star= np.array(f['a'])
with h5py.File('cen_po_star_gx.h5') as f:
    cen_po_star_gx = np.array(f['a'])    
with h5py.File('cen_po_gas_music.h5') as f:
    cen_po_gas= np.array(f['a'])
with h5py.File('cen_po_gas_gx.h5') as f:
    cen_po_gas_gx = np.array(f['a'])    
with h5py.File('GX_gas.h5') as f:
    dens_gas_gx = np.array(f['a'])
with h5py.File('GX_star.h5') as f:
    dens_star_gx = np.array(f['a'])
with h5py.File('music_gas.h5') as f:
    dens_gas = np.array(f['a'])
with h5py.File('music_star.h5') as f:
    dens_star = np.array(f['a'])
#### the 4th,figure the distribution out
plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'rainbow',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
plt.colorbar(label = r'$log[Number density]$')
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.scatter(cen_po_gas[0][0],cen_po_gas[0][1],s = 50,c = 'k',marker = 'x',)
plt.scatter(cen_po_gas[1][0],cen_po_gas[1][1],s = 50,c = 'k',marker = 'x',)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.savefig('MUSIC gas_cen_point No1_z_128',dpi=600)
plt.show()
plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
plt.colorbar(label = r'$log[Number density]$')
plt.scatter(po_halo[:,0],po_halo[:,1],s = 10,c = 'b',marker = 'o',alpha = 0.5)
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.savefig('MUSIC DM_cen_point No1_z_128',dpi=600)
plt.show()
if snap_N[-2] !=0:
    plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.scatter(cen_po_star[0][0],cen_po_star[0][1],s = 50,c = 'm',marker = 'o',)
    plt.scatter(cen_po_star[1][0],cen_po_star[1][1],s = 50,c = 'm',marker = 'o',)
    plt.xlim(4.98e5,5.028e5)
    plt.ylim(4.98e5,5.028e5)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC star_cen_point No1_z_128',dpi=600)
    plt.show()

plt.figure(figsize=(8,8))
plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
if snap_N[-2] !=0:
    plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.scatter(x0,y0,s = 50,c = 'k',marker = 'x',)
plt.scatter(po_halo[:,0],po_halo[:,1],s = 25,c = 'm',marker = 'o',alpha = 0.5)
plt.scatter(cen_po_star[0][0],cen_po_star[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star[1][0],cen_po_star[1][1],s = 50,c = 'm',marker = 'o',)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetMUSIC$')
plt.savefig('MUSIC Multi_point No1_z_128',dpi=600)
plt.show()

plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
if snap_N[-2] !=0:
    plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
plt.scatter(x0,y0,s = 50,c = 'k',marker = 'x',)
plt.scatter(cen_po_star[0][0],cen_po_star[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star[1][0],cen_po_star[1][1],s = 50,c = 'm',marker = 'o',)
plt.xlim(5.001e5,5.007e5)
plt.ylim(5.002e5,5.007e5)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetMUSIC$')
plt.savefig('MUSIC Multi_cen_comparation No1_z_128',dpi=600)
plt.show()
#############for GX
plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
           cmap = 'rainbow',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
plt.colorbar(label = r'$log[Number density]$')
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.scatter(cen_po_gas_gx[0][0],cen_po_gas_gx[0][1],s = 50,c = 'k',marker = 'x',)
plt.scatter(cen_po_gas_gx[1][0],cen_po_gas_gx[1][1],s = 50,c = 'k',marker = 'x',)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.savefig('GX gas_cen_point No1_z_128',dpi=600)
plt.show()
plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
plt.colorbar(label = r'$log[Number density]$')
plt.scatter(po_halo_gx[:,0],po_halo_gx[:,1],s = 10,c = 'b',marker = 'o',alpha = 0.5)
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.savefig('GX DM_cen_point No1_z_128',dpi=600)
plt.show()
if snap_N_gx[-2] !=0:
    plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.scatter(cen_po_star_gx[0][0],cen_po_star_gx[0][1],s = 50,c = 'm',marker = 'o',)
    plt.scatter(cen_po_star_gx[1][0],cen_po_star_gx[1][1],s = 50,c = 'm',marker = 'o',)
    plt.xlim(4.98e5,5.028e5)
    plt.ylim(4.98e5,5.028e5)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX star_cen_point No1_z_128',dpi=600)
    plt.show()

plt.figure(figsize=(8,8))
plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
if snap_N_gx[-2] !=0:
    plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.scatter(x0_gx,y0_gx,s = 50,c = 'k',marker = 'x',)
plt.scatter(po_halo_gx[:,0],po_halo_gx[:,1],s = 25,c = 'm',marker = 'o',alpha = 0.5)
plt.scatter(cen_po_star_gx[0][0],cen_po_star_gx[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star_gx[1][0],cen_po_star_gx[1][1],s = 50,c = 'm',marker = 'o',)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetX$')
plt.savefig('GX Multi_point No1_z_128',dpi=600)
plt.show()

plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
if snap_N_gx[-2] !=0:
    plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
plt.scatter(x0_gx,y0_gx,s = 50,c = 'k',marker = 'x',)
plt.scatter(cen_po_star_gx[0][0],cen_po_star_gx[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star_gx[1][0],cen_po_star_gx[1][1],s = 50,c = 'm',marker = 'o',)
plt.xlim(5.001e5,5.007e5)
plt.ylim(5.002e5,5.007e5)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetX$')
plt.savefig('GX Multi_cen_comparation No1_z_128',dpi=600)
plt.show()
#comparation
plt.figure(figsize = (8,8))
plt.subplot(1,2,1)
plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
if snap_N[-2] !=0:
    plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.scatter(x0,y0,s = 50,c = 'k',marker = 'x',)
plt.scatter(po_halo[:,0],po_halo[:,1],s = 25,c = 'm',marker = 'o',alpha = 0.5)
plt.scatter(cen_po_star[0][0],cen_po_star[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star[1][0],cen_po_star[1][1],s = 50,c = 'm',marker = 'o',)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetMUSIC$')
plt.subplot(1,2,2)
plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
if snap_N_gx[-2] !=0:
    plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
plt.xlim(4.98e5,5.028e5)
plt.ylim(4.98e5,5.028e5)
plt.scatter(x0_gx,y0_gx,s = 50,c = 'k',marker = 'x',)
plt.scatter(po_halo_gx[:,0],po_halo_gx[:,1],s = 25,c = 'm',marker = 'o',alpha = 0.5)
plt.scatter(cen_po_star_gx[0][0],cen_po_star_gx[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star_gx[1][0],cen_po_star_gx[1][1],s = 50,c = 'm',marker = 'o',)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetX$')
plt.tight_layout()
plt.savefig('Image halo_galaxy',dpi=600)
plt.show()
#central comparation
plt.figure(figsize = (8,8))
plt.subplot(1,2,1)
plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
if snap_N[-2] !=0:
    plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
plt.scatter(x0,y0,s = 50,c = 'k',marker = 'x',)
plt.scatter(cen_po_star[0][0],cen_po_star[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star[1][0],cen_po_star[1][1],s = 50,c = 'm',marker = 'o',)
plt.xlim(5.001e5,5.007e5)
plt.ylim(5.002e5,5.007e5)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetMUSIC$')
plt.subplot(1,2,2)
plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
if snap_N_gx[-2] !=0:
    plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
plt.scatter(x0_gx,y0_gx,s = 50,c = 'k',marker = 'x',)
plt.scatter(cen_po_star_gx[0][0],cen_po_star_gx[0][1],s = 50,c = 'm',marker = 'o',)
plt.scatter(cen_po_star_gx[1][0],cen_po_star_gx[1][1],s = 50,c = 'm',marker = 'o',)
plt.xlim(5.001e5,5.007e5)
plt.ylim(5.002e5,5.007e5)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetX$')
plt.tight_layout()
plt.savefig('Image central point',dpi=600)
plt.show()