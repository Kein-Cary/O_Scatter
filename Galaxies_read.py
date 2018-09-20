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
_id_ = 0
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
Mvir = np.array(main_halo['col4'])
Mstar = np.array(main_halo['col65'])
x0 = xhalo[0]
y0 = yhalo[0]
z0 = zhalo[0]
R0 = Rvir[0]
#try to find the subhalo belong to the first halo
halo_id = np.array(main_halo['col1'])
halo_id = np.int64(halo_id)
host_id = np.array(main_halo['col2'])
host_id = np.int64(host_id)
goal_halo = halo_id[_id_]
ix = host_id == goal_halo
po_halo = np.zeros((len(halo_id),3),dtype = np.float)
sub_mass = np.zeros(len(halo_id),dtype = np.float)
sub_star = np.zeros(len(halo_id),dtype = np.float)
for k in range(len(ix)):
    if k == _id_ :
        po_halo[k,:] = np.array([xhalo[k],yhalo[k],zhalo[k]])
        sub_mass[k] = Mvir[k]
        sub_star[k] = Mstar[k]
    elif ix[k] == True :
        po_halo[k,:] = np.array([xhalo[k],yhalo[k],zhalo[k]])
        sub_mass[k] = Mvir[k]
        sub_star[k] = Mstar[k]
    else:
        po_halo[k,:] = np.array([np.inf,np.inf,np.inf])
        sub_mass[k] = np.inf
        sub_star[k] = np.inf
sub_star = sub_star[sub_star != np.inf]
sub_mass = sub_mass[sub_mass != np.inf]
po_halo = po_halo[po_halo != np.inf]
po_halo = po_halo.reshape(np.int(len(po_halo)/3),3)
#isub = sub_star != 0
#isub = sub_star >= 6.5*10**10
isub = sub_mass >= 4.5*10**11
submass = sub_mass[isub]
substar = sub_star[isub]
pohalo = po_halo[isub,:]
#the above 7 lines is try to point the satillite galaxies with subhalo center
###next : try to point satillite galaxies with the center of star particles' density center
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

plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
if snap_N[-2] !=0:
    plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
plt.xlim(x0-R0,x0+R0)
plt.ylim(z0-R0,z0+R0)
plt.scatter(x0,z0,s = 50,c = 'k',marker = '+',)#point the center of the halo
plt.scatter(pohalo[:,0],pohalo[:,1],s = 20,c = 'm',marker = 'o',alpha = 0.5)
#point the satillite galaxies with the center of subhalo
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetMUSIC$')
plt.savefig('MUSIC galaxy_point in subhalo',dpi=600)
plt.show()
# next : try to find the center of number density
hist_s,edge_star = np.histogramdd(inl_star, bins = (15,15,15))
bin_x_star = edge_star[0]
bin_y_star = edge_star[1]
bin_z_star = edge_star[2]
ia = hist_s >= 50 
hists = hist_s[ia]
binx = np.zeros((len(hists),2),dtype = np.float)
biny = np.zeros((len(hists),2),dtype = np.float)
binz = np.zeros((len(hists),2),dtype = np.float)
#selection : choose those more densitive part to the next zoom in
inlstar = inl_star
hist_u = hist_s
import find 
for t in range(len(hists)):
    for p in range(hist_s.shape[0]):
        A = hist_u[p]
        B = hists[t] in A
        if B == True :
            ip_da = find.find2d(A,hists[t])
            ip_da = np.int0(ip_da)
            binx[t,:] = np.array([bin_x_star[p],bin_x_star[p+1]])
            biny[t,:] = np.array([bin_y_star[ip_da[0]],bin_y_star[ip_da[0]+1]])
            binz[t,:] = np.array([bin_z_star[ip_da[1]],bin_z_star[ip_da[1]+1]])
            hist_u[p,ip_da[0],ip_da[1]] = False
            break
        else:
            continue
cen_dens = np.zeros(len(hists),dtype = np.int0)
cen_po_star = np.zeros((len(hists),3),dtype = np.float)
for k in range(len(hists)):
    iv = ((inl_star[:,0]<=binx[k,1])&(inl_star[:,0]>=binx[k,0]))&((inl_star[:,1]<=biny[k,1])&(inl_star[:,1]>=biny[k,0])
            )&((inl_star[:,2]<=binz[k,1])&(inl_star[:,2]>=binz[k,0]))
    inlstar = inl_star[iv,:]
    hist_k,dege_k = np.histogramdd(inlstar, bins = (100,100,100))
    is_max = np.unravel_index(np.argmax(hist_k, axis=None), hist_k.shape)
    cen_dens[k] = hist_k[is_max]
    dege_k = np.array(dege_k)
    cen_x_star = (dege_k[0,is_max[0]+1]+dege_k[0,is_max[0]])/2.0
    cen_y_star = (dege_k[1,is_max[1]+1]+dege_k[1,is_max[1]])/2.0
    cen_z_star = (dege_k[2,is_max[2]+1]+dege_k[2,is_max[2]])/2.0
    cen_po_star[k,:] = np.array([cen_x_star,cen_y_star,cen_z_star])
#figure the results
plt.hist2d(inl_DM[:,0],inl_DM[:,1],bins = [1000,1000],
           cmap = 'viridis',vmin = 1e-1,vmax = snap_N[1]/100,norm=mpl.colors.LogNorm())
plt.hist2d(inl_gas[:,0],inl_gas[:,1],bins = [1000,1000],
           cmap = 'cool',vmin = 1e-1,vmax = snap_N[0]/100,norm=mpl.colors.LogNorm())
if snap_N[-2] !=0:
    plt.hist2d(inl_star[:,0],inl_star[:,1],bins = [1000,1000],
           cmap = 'plasma',vmin = 1e-1,vmax = snap_N[-2]/100,norm=mpl.colors.LogNorm())
plt.xlim(x0-R0,x0+R0)
plt.ylim(z0-R0,z0+R0)
plt.scatter(x0,z0,s = 50,c = 'k',marker = 'x',)#point the center of the halo
plt.scatter(cen_po_star[:,0],cen_po_star[:,1],s = 25,c = 'm',marker = 'o',alpha = 0.5)
#point the satillite galaxies with the centre of density(using starparticle)
plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.title(r'$GadgetMUSIC$')
plt.savefig('MUSIC galaxy_point in star',dpi=600)
plt.show()
