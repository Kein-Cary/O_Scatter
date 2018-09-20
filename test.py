# this file try to find the BCG in a given halo
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
dgas = np.sqrt((snap_shot_gas[:, 0]-x0)**2 +
               (snap_shot_gas[:, 1]-y0)**2 + (snap_shot_gas[:, 2]-z0)**2)
ig = dgas <= R0
inlgas = snap_shot_gas[ig, :]
inlmass_gas = snap_mass_gas[ig]

dDM = np.sqrt((snap_shot_DM[:, 0]-x0)**2 +
              (snap_shot_DM[:, 1]-y0)**2 + (snap_shot_DM[:, 2]-z0)**2)
iD = dDM <= R0
inlDM = snap_shot_DM[iD, :]
ddisk = np.sqrt((snap_shot_disk[:, 0]-x0)**2 +
                (snap_shot_disk[:, 1]-y0)**2 + (snap_shot_disk[:, 2]-z0)**2)
idisk = ddisk <= R0
inldisk = snap_shot_disk[idisk, :]
dbulge = np.sqrt((snap_shot_bulge[:, 0]-x0)**2 +
                 (snap_shot_bulge[:, 1]-y0)**2 + (snap_shot_bulge[:, 2]-z0)**2)
ibu = dbulge <= R0
inlbulge = snap_shot_bulge[ibu, :]

if snap_N[-2] != 0:
    dstar = np.sqrt((snap_shot_star[:, 0]-x0)**2+(snap_shot_star[:, 1]-y0)**2 + (snap_shot_star[:, 2]-z0)**2)
    ids = dstar <= R0
    inlstar = snap_shot_star[ids, :]
    inlmass_star = snap_mass_star[ids]
if snap_N[-1] != 0:
    dbndry = np.sqrt(
        (snap_shot_bndry[:, 0]-x0)**2+(snap_shot_bndry[:, 1]-y0)**2 + (snap_shot_bndry[:, 2]-z0)**2)
    idb = dbndry <= R0
    inlbndry = snap_shot_bndry[idb, :]
    inlmass_BH = snap_mass_BH[idb]
R_range = 50.0
edge_x = np.array([x0-R_range,x0+R_range])
edge_y = np.array([y0-R_range,y0+R_range])
edge_z = np.array([z0-R_range,z0+R_range])
iv_s = ((inlstar[:,0]<=edge_x[1])&(inlstar[:,0]>=edge_x[0]))&(
        (inlstar[:,1]<=edge_y[1])&(inlstar[:,1]>=edge_y[0]))&((inlstar[:,2]<=edge_z[1])&(inlstar[:,2]>=edge_z[0]))
test_star = inlstar[iv_s,:]
### to find the BCG with the way of KDtree
from sklearn.neighbors import KDTree as snKD
