import matplotlib as mpl
mpl.use('Agg')
import glob
import matplotlib.pyplot as plt

import h5py
import numpy as np
import astropy.units as U
import astropy.constants as C
from astropy import cosmology as apcy
from scipy.optimize import leastsq, curve_fit

## cosmology model
vc = C.c.to(U.km/U.s).value
Test_model = apcy.Planck15.clone(H0 = 67.74, Om0 = 0.311)
H0 = Test_model.H0.value
h = H0/100
Omega_m = Test_model.Om0
Omega_lambda = 1.-Omega_m
Omega_k = 1.- (Omega_lambda + Omega_m)
DH = vc/H0
G = C.G.value
Ms = C.M_sun.value
Mpc2km = U.Mpc.to(U.km)
Ms2kg = U.M_sun.to(U.kg)

## cluster catalogue
C_id = []
load1 = '/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/'
file1 = glob.glob('/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_*')
C_id = [str(file1[-4:]) for file1 in file1]
CID = np.sort(C_id)

def rho_t0(r, Mc, c):
	R = r
	c = c
	M = 10**Mc
	rho0 = (Mpc2km/Ms2kg)*(3*H0**2)/(8*np.pi*G)
	r200 = (3*M/(4*np.pi*rho0*200))**(1/3) 
	rs = r200/c
	rho0c = M/((np.log(1+c)-c/(1+c))*4*np.pi*rs**3)
	f0 = (R/rs)*(1+R/rs)**2
	f2 = rho0c/f0
	return f2, r200, rs

def rho_t(r, Mc, c):
	R = r
	c = c
	M = 10**Mc
	rho0 = (Mpc2km/Ms2kg)*(3*H0**2)/(8*np.pi*G)
	r200 = (3*M/(4*np.pi*rho0*200))**(1/3) 
	rs = r200/c
	rho0c = M/((np.log(1+c)-c/(1+c))*4*np.pi*rs**3)
	f0 = (R/rs)*(1+R/rs)**2
	f1 = rho0c/f0
	return f1

def d_rho(p, data, x):
	A, B = p
	y = rho_t(x, A, B)
	return y - data

def rho_fit():
	load2 = 'NewMDCLUSTER_'

	#sub_set = ['0120', '0135', 0153']
	## with params range: ([13.95, 1], [14.95, 9.95])
	#sub_set = ['0080', '0099', '0118', '0119'] ## with params range: ([13.95, 1], [14.95, 9.95])
	#sub_set = ['0132'] ## with params range: ([13.95, 1], [14.95, 9.95])
	#sub_set = ['0146'] ## with params range: ([13.95, 1], [15.05, 9.95])
	sub_set = ['0106', '0110', '0111', '0121', '0145', '0182']

	Nh = len(sub_set)
	for k in range(Nh):
		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/Cluster_profile/halo_tot_profile_%s.h5' % sub_set[k]) as f:
			halo_tot = np.array(f['a'])
		Rh = halo_tot[3,:]
		rho_g = halo_tot[2,:]
		rho_s = halo_tot[1,:]
		rho_d = halo_tot[0,:]

		rh = Rh[0,:][rho_d[0,:] != 0]
		rhoD = rho_d[0,:][rho_d[0,:] != 0]
		low_lim = 1e2
		rhoD = rhoD[rh >= low_lim]
		rh = rh[rh >= low_lim]
		'''
		Mc0 = 14
		C0 = 5
		fit_1 = leastsq(d_rho, [Mc0, C0], args = (rhoD[1:], rh[1:]))[0]
		mc = fit_1[0]
		cc = fit_1[1]
		'''
		popt, pcov = curve_fit(rho_t, rh, rhoD, bounds = ([13.95, 1], [14.95, 9.95]))
		mc = popt[0]
		cc = popt[1]

		print(mc)
		print(cc)
		fit_line, R_200, R_s = rho_t0(rh, mc, cc)

		plt.figure(figsize = (16, 9))
		plt.title('NFW fit for Cluster %s' % sub_set[k])
		plt.plot(rh, rhoD, 'b*', label = r'$Gadget-X \; \rho_{DM}$', alpha = 0.5)
		plt.plot(rh, fit_line, 'r--', label = r'$NFW \; fit$', alpha = 0.5)
		plt.axvline(x = R_200, color = 'k', ls = '-', label = r'$R_{200}$', alpha = 0.5)
		plt.axvline(x = R_s, color = 'k', ls = '--', label = r'$R_{s}$', alpha = 0.5)
		plt.xlim(low_lim, 2*np.max(rh))
		plt.xscale('log')
		plt.xlabel(r'$R[kpc / h]$')
		plt.yscale('log')
		plt.text(1.5 * low_lim, 1e5, s = r'$M_{h} = %.2fM_{\odot}/h$' % mc + '\n' + r'$C = %.2f$' % cc, fontsize = 'x-large')
		plt.ylabel(r'$\rho [h^{2} M_{\odot} / kpc^{3}]$')
		plt.legend(loc = 3, fontsize = 12)
		plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/compare_fig/NFW_fit_%s.png' 
			% sub_set[k], dpi = 300)
		plt.close()
	raise
	return

def main():
	rho_fit()

if __name__ == "__main__":
	main()
