#this file try to read the snapshot data to find the particle information
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
#section1: read out the scale factor and redshift
f1 = asc.read('D:/mask/snapshot/MUSIC/redshift_list.txt')
scale_factor = f1['a']
redshift = f1['z']
#next try to sort redshift and scale factor and read snap_shot as redshift order
a = np.array(scale_factor)
a = a[::-1]
z = np.array(redshift)
z = z[::-1]
#set the goal snap_shot
No_snap = '128'
_id_ = 0
#section2:read out the position data of MUSIC
snap_z = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'redshift')
snap_z = np.abs(snap_z)
print('Now redshift is %.3f'%snap_z)
snap_N = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'npartTotal')
print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
      %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
snap_name = pygdr.readheader('D:/mask/snapshot/MUSIC/snap_%s'%No_snap,'header')
raise
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
#the total mass distribution of snap_shot_128
##next,try to find the particles belong to halo 128000000000001(in MUSIC simulation)
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
#section3:read out the position data of GadgetX
#########for GX simulation
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
x0_gx = xhalo_gx[_id_]
y0_gx = yhalo_gx[_id_]
z0_gx = zhalo_gx[_id_]
R0_gx = Rvir_gx[_id_]
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

def fig_music(t):
    '''
    if snap_N[-2] !=0:
        plt.hist2d(snap_shot_star[:,0],snap_shot_star[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('MUSIC star distribution z_128',dpi=600)
        plt.show()
    plt.hist2d(snap_shot_gas[:,0],snap_shot_gas[:,1],bins = [1000,1000],
               cmap = 'rainbow',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC gas distribution z_128',dpi=600)
    plt.show()
    plt.hist2d(snap_shot_DM[:,0],snap_shot_DM[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC DM distribution z_128',dpi=600)
    plt.show()
    '''
    '''
    if snap_N[-1] !=0:
        plt.hist2d(snap_shot_bndry[:,0],snap_shot_bndry[:,1],bins = [500,500],
                   cmap = 'GnBu',vmin = 1e-1,vmax = snap_N[-1]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('MUSIC boundary distribution z_128',dpi=600)
        plt.show()
    plt.hist2d(snap_shot_disk[:,0],snap_shot_disk[:,1],bins = [500,500],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N[2]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC disk distribution z_128',dpi=600)
    plt.show()
    plt.hist2d(snap_shot_bulge[:,0],snap_shot_bulge[:,1],bins = [500,500],
               cmap = 'summer',vmin = 1e-1,vmax = snap_N[3]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC bulge distribution z_128',dpi=600)
    plt.show()
    '''
    #for main halo
    plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'rainbow',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.xlim(x0-R0,x0+R0)
    plt.ylim(y0-R0,y0+R0)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC gas distribution No1_z_128',dpi=600)
    plt.show()

    plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0-R0,x0+R0)
    plt.ylim(y0-R0,y0+R0)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC DM distribution No1_z_128',dpi=600)
    plt.show()
    if snap_N[-2] !=0:
        plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        #plt.xlim(4.98e5,5.028e5)
        #plt.ylim(4.98e5,5.028e5)
        plt.xlim(x0-R0,x0+R0)
        plt.ylim(y0-R0,y0+R0)
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('MUSIC star distribution No1_z_128',dpi=600)
        plt.show()
    '''
    plt.hist2d(inl_disk[:,0],inl_disk[:,1],bins = [500,500],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N[2]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0-R0,x0+R0)
    plt.ylim(y0-R0,y0+R0)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC disk distribution No1_z_128',dpi=600)
    plt.show()
    plt.hist2d(inl_bulge[:,0],inl_bulge[:,1],bins = [500,500],
               cmap = 'summer',vmin = 1e-1,vmax = snap_N[3]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0-R0,x0+R0)
    plt.ylim(y0-R0,y0+R0)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('MUSIC bulge distribution No1_z_128',dpi=600)
    plt.show()
    if snap_N[-1] !=0:
        plt.hist2d(inl_bndry[:,0],inl_bndry[:,1],bins = [500,500],
               cmap = 'gist_rainbow',vmin = 1e-1,vmax = snap_N[-1]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        #plt.xlim(4.98e5,5.028e5)
        #plt.ylim(4.98e5,5.028e5)
        plt.xlim(x0-R0,x0+R0)
        plt.ylim(y0-R0,y0+R0)
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('MUSIC bndry distribution No1_z_128',dpi=600)
        plt.show()
    '''
    #multi-distribution
    plt.figure(figsize=(8,8))
    plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0-R0,x0+R0)
    plt.ylim(y0-R0,y0+R0)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetMUSIC$')
    plt.savefig('MUSIC Multi-distribution No1_z_128',dpi=600)
    plt.show()
    plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
    plt.scatter(x0,y0,s = 50,c = 'k',marker = 'x',)
    #plt.xlim(5.001e5,5.007e5)
    #plt.ylim(5.002e5,5.007e5)
    plt.xlim(x0-R0*3.5/22.0,x0+R0*3.5/22.0)
    plt.ylim(y0-R0*3.0/14.2,y0+R0*3.0/14.2)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetMUSIC$')
    plt.savefig('MUSIC Multi_center No1_z_128',dpi=600)
    plt.show()
    return
fig_music(t = True)
def fig_GX(p):
    '''
    if snap_N_gx[-2] !=0:
        plt.hist2d(snap_shot_star_gx[:,0],snap_shot_star_gx[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('GX star distribution z_128',dpi=600)
        plt.show()
    plt.hist2d(snap_shot_gas_gx[:,0],snap_shot_gas_gx[:,1],bins = [1000,1000],
               cmap = 'rainbow',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX gas distribution z_128',dpi=600)
    plt.show()
    plt.hist2d(snap_shot_DM_gx[:,0],snap_shot_DM_gx[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX DM distribution z_128',dpi=600)
    plt.show()
    '''
    '''
    if snap_N_gx[-1] !=0:
        plt.hist2d(snap_shot_bndry_gx[:,0],snap_shot_bndry_gx[:,1],bins = [500,500],
                   cmap = 'GnBu',vmin = 1e-1,vmax = snap_N_gx[-1]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('GX boundary distribution z_128',dpi=600)
        plt.show()
    plt.hist2d(snap_shot_disk_gx[:,0],snap_shot_disk_gx[:,1],bins = [500,500],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[2]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX disk distribution z_128',dpi=600)
    plt.show()
    plt.hist2d(snap_shot_bulge_gx[:,0],snap_shot_bulge_gx[:,1],bins = [500,500],
               cmap = 'summer',vmin = 1e-1,vmax = snap_N_gx[3]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX bulge distribution z_128',dpi=600)
    plt.show()
    '''
    #for main halo
    plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
               cmap = 'rainbow',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
    plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX gas distribution No1_z_128',dpi=600)
    plt.show()
    plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
    plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX DM distribution No1_z_128',dpi=600)
    plt.show()
    if snap_N[-2] !=0:
        plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        #plt.xlim(4.98e5,5.028e5)
        #plt.ylim(4.98e5,5.028e5)
        plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
        plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('GX star distribution No1_z_128',dpi=600)
        plt.show()
    '''
    plt.hist2d(inl_disk_gx[:,0],inl_disk_gx[:,1],bins = [500,500],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[2]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
    plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX disk distribution No1_z_128',dpi=600)
    plt.show()
    plt.hist2d(inl_bulge_gx[:,0],inl_bulge_gx[:,1],bins = [500,500],
               cmap = 'summer',vmin = 1e-1,vmax = snap_N_gx[3]/100,norm=mpl.colors.LogNorm())
    plt.colorbar(label = r'$log[Number density]$')
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
    plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('GX bulge distribution No1_z_128',dpi=600)
    plt.show()
    if snap_N[-1] !=0:
        plt.hist2d(inl_bndry_gx[:,0],inl_bndry_gx[:,1],bins = [500,500],
               cmap = 'gist_rainbow',vmin = 1e-1,vmax = snap_N_gx[-1]/100,norm=mpl.colors.LogNorm())
        plt.colorbar(label = r'$log[Number density]$')
        #plt.xlim(4.98e5,5.028e5)
        #plt.ylim(4.98e5,5.028e5)
        plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
        plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.savefig('GX bndry distribution No1_z_128',dpi=600)
        plt.show()
    '''
    #multi-distribution
    plt.figure(figsize=(8,8))
    plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
    plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetX$')
    plt.savefig('GX Multi-distribution No1_z_128',dpi=600)
    plt.show()
    #central area view
    plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
    plt.scatter(x0_gx,y0_gx,s = 50,c = 'k',marker = 'x',)
    #plt.xlim(5.001e5,5.007e5)
    #plt.ylim(5.002e5,5.007e5)
    plt.xlim(x0_gx-R0_gx*3.5/22.0,x0_gx+R0_gx*3.5/22.0)
    plt.ylim(y0_gx-R0_gx*3.0/14.2,y0_gx+R0_gx*3.0/14.2)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetX$')
    plt.savefig('GX Multi-center No1_z_128',dpi=600)
    plt.show()
    return
fig_GX(p = True)
def fig_comparation(q):
    #comparation from x-z panel
    plt.figure(figsize = (8,8))
    plt.subplot(1,2,1)
    plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0-R0,x0+R0)
    plt.ylim(y0-R0,y0+R0)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetMUSIC$')
    plt.subplot(1,2,2)
    plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
    #plt.xlim(4.98e5,5.028e5)
    #plt.ylim(4.98e5,5.028e5)
    plt.xlim(x0_gx-R0_gx,x0_gx+R0_gx)
    plt.ylim(y0_gx-R0_gx,y0_gx+R0_gx)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetX$')
    plt.tight_layout()
    plt.savefig('Image halo particle',dpi=600)
    plt.show()
    #for central area
    plt.figure()
    plt.subplot(1,2,1)
    plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
    plt.scatter(x0,y0,s = 50,c = 'k',marker = 'x',)
    #plt.xlim(5.001e5,5.007e5)
    #plt.ylim(5.002e5,5.007e5)
    plt.xlim(x0-R0*3.5/22.0,x0+R0*3.5/22.0)
    plt.ylim(y0-R0*3.0/14.2,y0+R0*3.0/14.2)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetMUSIC$')
    plt.subplot(1,2,2)
    plt.hist2d(inl_DM_gx[:,0],inl_DM_gx[:,1],bins = [1000,1000],
               cmap = 'viridis',vmin = 1e-1,vmax = snap_N_gx[1]/100,norm=mpl.colors.LogNorm())
    plt.hist2d(inl_gas_gx[:,0],inl_gas_gx[:,1],bins = [1000,1000],
               cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[0]/100,norm=mpl.colors.LogNorm())
    if snap_N[-2] !=0:
        plt.hist2d(inl_star_gx[:,0],inl_star_gx[:,1],bins = [1000,1000],
               cmap = 'plasma',vmin = 1e-1,vmax = snap_N_gx[-2]/100,norm=mpl.colors.LogNorm())
    plt.scatter(x0_gx,y0_gx,s = 50,c = 'k',marker = 'x',)
    #plt.xlim(5.001e5,5.007e5)
    #plt.ylim(5.002e5,5.007e5)
    plt.xlim(x0_gx-R0_gx*3.5/22.0,x0_gx+R0_gx*3.5/22.0)
    plt.ylim(y0_gx-R0_gx*3.0/14.2,y0_gx+R0_gx*3.0/14.2)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])
    plt.title(r'$GadgetX$')
    plt.tight_layout()
    plt.savefig('Image central',dpi=600)
    plt.show()
    return
fig_comparation(q = True)