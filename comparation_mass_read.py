#use the file data(total_mass*.h5)to figure
import matplotlib as mpl
import astropy.io.ascii as asc
import numpy as np
import h5py
import pandas as pd
from handy import scatter as hsc
import matplotlib.pyplot as plt
#'D:/mask/MUSIC/MUSIC_reshift/NewMDCLUSTER_0001/'
#'D:/mask/G_X/G_x_redshift/NewMDCLUSTER_0001/'
with h5py.File('total_mass_MUSIC.h5') as f:
    y0 = f['a']
    mass_MUSIC = np.array(y0)
with h5py.File('total_mass_GX.h5') as f:
    y1 = f['a']
    mass_GX = np.array(y1)
with h5py.File('D:/python1/pydocument/O_Scatter/MUSIC_reshift/NewMDCLUSTER_0001/redshift.h5') as f:
    y2 = f['a']
    redshift = np.array(y2)
with h5py.File('D:/python1/pydocument/O_Scatter/G_x_redshift/NewMDCLUSTER_0001/Redshift_GX.h5') as f:
    y3 = f['a']
    redshift_gx = np.array(y3)
#data of MUSIC simulation
halo_mass = np.array(mass_MUSIC[0,:])
star_mass = np.array(mass_MUSIC[1,:])
gas_mass = np.array(mass_MUSIC[2,:])
baryon_m = np.array(mass_MUSIC[3,:])
dark_mass = np.array(mass_MUSIC[4,:])
#data of GadgetX simulation
halo_mass_gx = np.array(mass_GX[0,:])
star_mass_gx = np.array(mass_GX[1,:])
gas_mass_gx = np.array(mass_GX[2,:])
baryon_m_gx = np.array(mass_GX[3,:])
dark_mass_gx = np.array(mass_GX[4,:])
#figure
ia = redshift<=1.0 
z_MUSIC = np.array(redshift[ia])
ia_gx = redshift_gx<=1.0 
z_GX = np.array(redshift_gx[ia_gx])
#find the mass change in a short time
N = len(z_MUSIC[0:-1])
M = len(z_GX[0:-1])
deta_mh = np.zeros(N,dtype = np.float)
deta_ms = np.zeros(N,dtype = np.float)
deta_mg = np.zeros(N,dtype = np.float)
deta_mb = np.zeros(N,dtype = np.float)
deta_md = np.zeros(N,dtype = np.float)
#for GadgetX 
deta_mh_gx = np.zeros(M,dtype = np.float)
deta_ms_gx = np.zeros(M,dtype = np.float)
deta_mg_gx = np.zeros(M,dtype = np.float)
deta_mb_gx = np.zeros(M,dtype = np.float)
deta_md_gx = np.zeros(M,dtype = np.float)
for k in range(N):
    dz = z_MUSIC[k+1]-z_MUSIC[k]
    deta_mh[k] = (halo_mass[k]-halo_mass[k+1])/dz
    deta_ms[k] = (star_mass[k]-star_mass[k+1])/dz
    deta_mg[k] = (gas_mass[k]-gas_mass[k+1])/dz
    deta_mb[k] = (baryon_m[k]-baryon_m[k+1])/dz
    deta_md[k] = (dark_mass[k]-dark_mass[k+1])/dz
for k in range(M):
    dz_gx = z_GX[k+1]-z_GX[k]
    deta_mh_gx[k] = (halo_mass_gx[k]-halo_mass_gx[k+1])/dz_gx
    deta_ms_gx[k] = (star_mass_gx[k]-star_mass_gx[k+1])/dz_gx
    deta_mg_gx[k] = (gas_mass_gx[k]-gas_mass_gx[k+1])/dz_gx
    deta_mb_gx[k] = (baryon_m_gx[k]-baryon_m_gx[k+1])/dz_gx
    deta_md_gx[k] = (dark_mass_gx[k]-dark_mass_gx[k+1])/dz_gx
plt.figure(figsize=(8,9))
plt.plot(z_MUSIC[0:-1],deta_ms,'r-',label = r'$dM_\ast/dz$')
plt.plot(z_MUSIC[0:-1],deta_mg,'g-',label = r'$dM_g/dz$')
plt.plot(z_MUSIC[0:-1],deta_mb,'m-',label = r'$dM_b/dz$')
plt.plot(z_MUSIC[0:-1],deta_md,'k-',label = r'$dM_d/dz$')

plt.plot(z_GX[0:-1],deta_ms_gx,'r--',label = r'$dM_\ast/dz -GX$')
plt.plot(z_GX[0:-1],deta_mg_gx,'g--',label = r'$dM_g/dz -GX$')
plt.plot(z_GX[0:-1],deta_mb_gx,'m--',label = r'$dM_b/dz -GX$')
plt.plot(z_GX[0:-1],deta_md_gx,'k--',label = r'$dM_d/dz -GX$')
plt.yscale('log')
#plt.legend(loc = 1)
plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=1,
   ncol=2, fontsize=12., borderaxespad=1.) 
plt.xlabel(r'$z$')
plt.ylabel(r'$M[M_\odot /h]$')
plt.savefig('mass growth',dpi = 600)
plt.show()

#set the goal range of redshift
plt.figure(figsize=(8,9))
plt.plot(z_MUSIC,halo_mass,'b-',label = r'$M_h[MUSIC]$')
plt.plot(z_MUSIC,star_mass,'r-',label = r'$M_s[MUSIC]$')
plt.plot(z_MUSIC,gas_mass,'g-',label = r'$M_g[MUSIC]$')
plt.plot(z_MUSIC,baryon_m,'m-',label = r'$M_b[MUSIC]$')
plt.plot(z_MUSIC,dark_mass,'k-',label = r'$M_d[MUSIC]$')
plt.plot(z_GX,halo_mass_gx,'b--',label = r'$M_h[GX]$')
plt.plot(z_GX,star_mass_gx,'r--',label = r'$M_s[GX]$')
plt.plot(z_GX,gas_mass_gx,'g--',label = r'$M_g[GX]$')
plt.plot(z_GX,baryon_m_gx,'m--',label = r'$M_b[GX]$')
plt.plot(z_GX,dark_mass_gx,'k--',label = r'$M_d[GX]$')
plt.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=2,
   ncol=5, mode = "expand",fontsize=10., borderaxespad=1.) 
plt.xlabel(r'$ z $')
plt.ylabel(r'$M [M_\odot /h]$')
plt.yscale('log')
plt.title('M-z relation')
plt.savefig('M-z',dpi=600)
plt.show()

f = plt.figure(figsize = (8,9))
f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
ratio_gas_b = gas_mass/baryon_m
ratio_gas_d = gas_mass/dark_mass
ratio_star_d = star_mass/dark_mass
ratio_star_b = star_mass/baryon_m
ratio_b = baryon_m/halo_mass
ratio_d = dark_mass/halo_mass
ratio_gas_b_gx = gas_mass_gx/baryon_m_gx
ratio_gas_d_gx = gas_mass_gx/dark_mass_gx
ratio_star_d_gx = star_mass_gx/dark_mass_gx
ratio_star_b_gx = star_mass_gx/baryon_m_gx
ratio_b_gx = baryon_m_gx/halo_mass_gx
ratio_d_gx = dark_mass_gx/halo_mass_gx

ax1.plot(z_MUSIC,ratio_gas_b,'b-',label = r'$\eta -M_g-M_b$')
ax1.plot(z_MUSIC,ratio_d,'k-',label = r'$\eta -M_d-M_h$')
ax1.plot(z_GX,ratio_gas_b_gx,'b--',label = r'$\eta -M_g-M_b -GX$')
ax1.plot(z_GX,ratio_d_gx,'k--',label = r'$\eta -M_d-M_h -GX$')
ax1.legend(bbox_to_anchor=(0., 0.85, 1., 0.), loc=2,
     ncol=2, fontsize=7.5, borderaxespad=0.) 
plt.ylabel(r'$\eta$')
plt.ylim(0.75,0.95)
ax1.set_title(r'$\eta - z $')

ax2.plot(z_MUSIC,ratio_gas_d,'g-',label = r'$\eta[M_g-M_d]$')
ax2.plot(z_MUSIC,ratio_star_d,'r-',label = r'$\eta[M_\ast-M_d]$')
ax2.plot(z_MUSIC,ratio_star_b,'m-',label = r'$\eta[M_\ast-M_b]$')
ax2.plot(z_MUSIC,ratio_b,'y-',label = r'$\eta[M_b-M_h]$')
ax2.plot(z_GX,ratio_gas_d_gx,'g--',label = r'$\eta[M_g-M_d]GX$')
ax2.plot(z_GX,ratio_star_d_gx,'r--',label = r'$\eta[M_\ast-M_d]GX$')
ax2.plot(z_GX,ratio_star_b_gx,'m--',label = r'$\eta[M_\ast-M_b]GX$')
ax2.plot(z_GX,ratio_b_gx,'y--',label = r'$\eta[M_b-M_h]GX$')
ax2.legend(bbox_to_anchor=(0., 1., 1., 0.), loc=2,
     ncol=4, fontsize=7.5, borderaxespad=0.) 
plt.ylim(0.,0.3)
plt.ylabel(r'$\eta$')
plt.xlabel(r'$ z $')
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.savefig('ratio of mass',dpi=600)
plt.show()
