import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob
import h5py
import numpy as np
import pandas as pd
import astropy.io.ascii as asc
import scipy.integrate as interg
import astropy.units as U
import astropy.constants as C

import find
import changds
import pygadgetreader as pygdr

## constant from Plank2015
h = 0.678
Omk = 0
Omm = 0.308
Mpc2km = U.Mpc.to(U.km)
H0 = 100*h

resolution = 1. # in unit kpc/h
BCG_size = 50 # in unit kpc/h
R_in = 50
R_out = 100

C_id = [] # cluster ID
load1 = '/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/'
file1 = glob.glob('/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_*')
C_id = [str(file1[-4:]) for file1 in file1]
CID = np.sort(C_id)

def rho_GX():
	load2 = 'NewMDCLUSTER_'
	sub_set = ['0145', '0153', '0182']
	'''
	sub_set = C_id
	'''
	Nh = len(sub_set)
	bins = 50
	for k in range(Nh):

		if (sub_set[k] == '0080') | (sub_set[k] == '0110'):
			trc_id = 1
		else:
			trc_id = 0

		sub_file = glob.glob(
			'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/*.AHF_halos' % sub_set[k])
		snap_id = [str(sub_file[107:110]) for sub_file in sub_file]
		snap = np.sort(snap_id)

		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/z_for_C_%s.h5' % sub_set[k]) as f:
			zc = np.array(f['a'])
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/cluster_%s_tree.h5' % sub_set[k]) as f:
			HID = np.array(f['a'])

		L = HID.shape[0]
		iv = np.zeros(L, dtype = np.int0)
		for ll in range(L):
			u = find.find1d(HID[ll,:],0)
			iv[ll] = u
		id_z = zc[iv[trc_id]]
		N = len(zc) - 1

		rho_s = np.zeros((N, 3*bins), dtype = np.float)
		rho_g = np.zeros((N, 3*bins), dtype = np.float)
		rho_dm = np.zeros((N, 3*bins), dtype = np.float)
		R_h = np.zeros((N, 3*bins), dtype = np.float)

		rho_bcg = np.zeros(bins, dtype = np.float)
		R_bcg = np.zeros(bins, dtype = np.float)
		# BCG + ICM component
		rho_ccm = np.zeros((N, 3*bins), dtype = np.float)
		R_ccm = np.zeros((N, 3*bins), dtype = np.float)

		for qq in range(2):
			zq = qq
			z_0 = zc[zq]
			Goal_ID = snap[N - zq]

			halo_list = asc.read(
				'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/GadgetX-NewMDCLUSTER_%s.snap_%s.z%.3f.AHF_halos'
				% (sub_set[k], sub_set[k], Goal_ID, z_0), converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
			Rvir = np.array(halo_list['col12'])
			xhalo = np.array(halo_list['col6'])
			yhalo = np.array(halo_list['col7'])
			zhalo = np.array(halo_list['col8'])
			Nbins = np.array(halo_list['col37'])
			Mhalo = np.array(halo_list['col4'])
			Mgas = np.array(halo_list['col45'])
			Mstar= np.array(halo_list['col65'])
			halo_colum = np.array(halo_list['col1'])
			host_ID = np.array(halo_list['col2'])

			goal_halo = HID[trc_id, qq]
			check_id = find.find1d(halo_colum, goal_halo)
			x0 = xhalo[check_id]
			y0 = yhalo[check_id]
			z0 = zhalo[check_id]
			R0 = Rvir[check_id]

			snap_N = pygdr.readheader(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'npartTotal')
			snap_shot_gas = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'pos', 'gas')
			snap_mass_gas = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'mass', 'gas')
			snap_shot_dm = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'pos', 'dm')
			snap_mass_dm = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'mass', 'dm')	
			try:
				snap_shot_star = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'pos', 'star')
				snap_mass_star = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'mass','star')
			except SystemExit:
				print('no star particles now')
				break
			# for ICM, BCG component
			try:
				dstar = np.sqrt((snap_shot_star[:, 0] - x0)**2 + 
					(snap_shot_star[:, 1] - y0)**2 + (snap_shot_star[:, 2] - z0)**2)
				ids = dstar <= R0
				inlstar = snap_shot_star[ids, :]
				inlmass_star = snap_mass_star[ids]

				R_range = 150.0
				edge_x = np.array([x0 - R_range, x0 + R_range])
				edge_y = np.array([y0 - R_range, y0 + R_range])
				edge_z = np.array([z0 - R_range, z0 + R_range])
				iv_s = ((inlstar[:,0] <= edge_x[1])&(inlstar[:,0] >= edge_x[0]))&(
						(inlstar[:,1] <= edge_y[1])&(inlstar[:,1] >= edge_y[0]))&(
						(inlstar[:,2] <= edge_z[1])&(inlstar[:,2] >= edge_z[0]))
				inl_star = inlstar[iv_s,:]
				num_bins = np.ceil(R_range*2 /resolution)
				try:
					hist_star, edge_star = np.histogramdd(inl_star, bins=(num_bins, num_bins, num_bins))
					bin_x_star = np.array(edge_star[0])
					bin_y_star = np.array(edge_star[1])
					bin_z_star = np.array(edge_star[2])

					inumber1 = hist_star >= 10
					maxN = len(inumber1)
					cen_po_star = np.zeros((maxN,3),dtype = np.float)
					hist_use = hist_star
					for p in range(maxN):
						is_max = np.unravel_index(np.argmax(hist_use, axis=None), hist_use.shape)
						cenxstar = (bin_x_star[is_max[0] + 1] + bin_x_star[is_max[0]])/2.0
						cenystar = (bin_y_star[is_max[1] + 1] + bin_y_star[is_max[1]])/2.0
						cenzstar = (bin_z_star[is_max[2] + 1] + bin_z_star[is_max[2]])/2.0
						cen_po_star[p,:] = np.array([cenxstar,cenystar,cenzstar])
						hist_use[is_max] = 0.0
					compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+
						(cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
					ismin = find.find1d(compare_d,np.min(compare_d))
					cen_x_star = cen_po_star[ismin, 0]
					cen_y_star = cen_po_star[ismin, 1]
					cen_z_star = cen_po_star[ismin, 2]
					## rho_bcg part
					R_BCG = BCG_size
					r_bcg = np.logspace(-3, np.log10(R_BCG), bins + 1)
					n_r = len(r_bcg)

					inl_m_s = np.zeros(n_r,dtype = np.float)
					r_star = np.sqrt((inlstar[:,0]-cen_x_star)**2+(inlstar[:,1]-cen_y_star)**2+
							(inlstar[:,2]-cen_z_star)**2)
					for p in range(n_r):
						ib = r_star <= r_bcg[p]
						inl_m_s[p] = np.sum(inlmass_star[ib]*10**10)
					for p in range(n_r-1):
						dr = r_bcg[p+1] - r_bcg[p]
						rho_bcg[p] = (inl_m_s[p+1]-inl_m_s[p]) / (4*np.pi*dr*r_bcg[p]**2)
						R_bcg[p] = 0.5 * (r_bcg[p+1] + r_bcg[p])					
					try:
						ih = host_ID == halo_colum[trc_id]
						msat = Mstar[ih]
						sub_x = xhalo[ih]
						sub_y = yhalo[ih]
						sub_z = zhalo[ih]
						sub_r = Rvir[ih]

						dr_sub0 = np.sqrt((sub_x - x0)**2 
							+ (sub_y - y0)**2 + (sub_z - z0)**2)
						isub = dr_sub0 <= R0
						real_sub = msat[isub]
						real_sub_x = sub_x[isub]
						real_sub_y = sub_y[isub]
						real_sub_z = sub_z[isub]
						real_sub_r = sub_r[isub]

						dr_sub1 = np.sqrt((real_sub_x - cen_x_star)**2 + 
							(real_sub_y - cen_y_star)**2 + (real_sub_z - cen_z_star)**2)
						ucen = np.where(dr_sub1 == np.min(dr_sub1))[0]
						cen_mass = real_sub[ucen[0]]
					except ValueError:
						print('there is noly one halo!')
					# select ICL + BCG from all star 
					drs = np.sqrt((inlstar[:, 0] - cen_x_star)**2 +
						(inlstar[:, 1] - cen_y_star)**2 +(inlstar[:, 2] - cen_z_star)**2)

					pos = inlstar * 1
					out_star = inlmass_star * 1
					r_icm = np.logspace(-3, np.log10(R0), 3*bins + 1)
					for tt in range(len(sub_r)):
						dr = np.sqrt((pos[:,0] - sub_x[tt])**2 + 
							(pos[:,1] - sub_y[tt])**2 + (pos[:,2] - sub_z[tt])**2)
						iic = dr >= sub_r[tt]
						pos = pos[iic, :]
						out_star = out_star[iic]
					dr = np.sqrt((pos[:,0] - x0)**2 + (pos[:,1] - y0)**2 +(pos[:,2] - z0)**2)
					dm_icm = np.zeros(len(r_icm), dtype = np.float)
					for kk in range(len(r_icm)):
						idr = dr <= r_icm[kk]
						ddm = out_star[idr]
						dm_icm[kk] = np.sum(ddm) * 10**10
					for kk in range(len(r_icm) - 1):
						ddr = r_icm[kk+1] - r_icm[kk]
						rho_ccm[qq, kk] = (dm_icm[kk + 1] - dm_icm[kk]) / (4*np.pi * ddr * r_icm[kk]**2)
						R_ccm[qq, kk] = 0.5 * (r_icm[kk + 1] + r_icm[kk])
				except ValueError:
					continue
			except IndexError:
				continue
			## halo profile
			dr_gas = np.sqrt((snap_shot_gas[:, 0] - x0)**2 + 
				(snap_shot_gas[:,1] - y0)**2 + (snap_shot_gas[:,2] - z0)**2)
			i_gas = dr_gas <= R0
			inl_gas = snap_shot_gas[i_gas,:]
			inl_m_gas = snap_mass_gas[i_gas]

			dr_star = np.sqrt((snap_shot_star[:, 0] - x0)**2 + 
				(snap_shot_star[:, 1] - y0)**2 + (snap_shot_star[:, 2] - z0)**2)
			i_star = dr_star <= R0
			inl_star = snap_shot_star[i_star,:]
			inl_m_star = snap_mass_star[i_star]

			dr_dm = np.sqrt((snap_shot_dm[:, 0] - x0)**2 + 
				(snap_shot_dm[:, 1] - y0)**2 + (snap_shot_dm[:, 2] - z0)**2)
			i_dm = dr_dm <= R0
			inl_dm = snap_shot_dm[i_dm, :]
			inl_m_dm = snap_mass_dm[i_dm]

			subr = np.logspace(-3, np.log10(R0), 3*bins + 1)
			sub_dr_d = np.sqrt((inl_dm[:, 0] - x0)**2 + 
				(inl_dm[:, 1] - y0)**2 + (inl_dm[:, 2] - z0)**2)
			sub_dr_s = np.sqrt((inl_star[:, 0] - x0)**2 + 
				(inl_star[:, 1] - y0)**2 + (inl_star[:, 2] - z0)**2)
			sub_dr_g = np.sqrt((inl_gas[:, 0] - x0)**2 + 
				(inl_gas[:, 1] - y0)**2 + (inl_gas[:, 2] - z0)**2)
			md = []
			ms = []
			mg = []
			for tt in range(len(subr)):
				idr = sub_dr_d < subr[tt]
				sub_dm = inl_m_dm[idr]
				subdm = np.sum(sub_dm) * 10**10
				md.append(subdm)

				idr = sub_dr_s < subr[tt]
				sub_star = inl_m_star[idr]
				substar = np.sum(sub_star) * 10**10
				ms.append(substar)

				idr = sub_dr_g <= subr[tt]
				sub_gas = inl_m_gas[idr]
				subgas = np.sum(sub_gas) * 10**10
				mg.append(subgas)

			md = np.array(md)
			ms = np.array(ms)
			mg = np.array(mg)
			for tt in range(len(subr) - 1):
				dr_d = subr[tt+1] - subr[tt]
				rho_s[tt] = (ms[tt + 1] - ms[tt]) / (4* np.pi* dr_d* subr[tt]**2)
				rho_g[tt] = (mg[tt + 1] - mg[tt]) / (4* np.pi* dr_d* subr[tt]**2)
				rho_dm[tt] = (md[tt + 1] - md[tt]) / (4* np.pi* dr_d* subr[tt]**2)
				mean_r = 0.5 * (subr[tt + 1] + subr[tt])
				R_h[tt] = mean_r * 1

		cord_array = np.array([rho_ccm, R_ccm])
		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/Cluster_profile/halo_bcg_icm_profile_%s.h5' % sub_set[k], 'w') as f:
			f['a'] = np.array(cord_array)
		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/Cluster_profile/halo_bcg_icm_profile_%s.h5' % sub_set[k] ) as f:
			for qq in range(len(cord_array)):
				f['a'][qq,:] = cord_array[qq,:]

		sum_array = np.array([rho_dm, rho_s, rho_g, R_h])
		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/Cluster_profile/halo_tot_profile_%s.h5' % sub_set[k], 'w') as f:
			f['a'] = np.array(sum_array)
		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/Cluster_profile/halo_tot_profile_%s.h5' % sub_set[k] ) as f:
			for tt in range(len(sum_array)):
				f['a'][tt, :] = sum_array[tt, :]

	return

def fig_out():
	load2 = 'NewMDCLUSTER_'
	sub_set = ['0153', '0182']
	'''
	sub_set = CID
	'''
	Nh = len(sub_set)
	bins = 50

	rh = []
	rhog = []
	rhos = []
	rhod = []

	rccm = []
	rhoccm = []

	etaccm = []
	for k in range(Nh):
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/z_for_C_%s.h5' % sub_set[k]) as f:
			zc = np.array(f['a'])
		zq = 0	
		z0 = zc[zq]

		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/halo_tot_profile_%s_z_%.3f.h5' % (sub_set[k], z0)) as f:
			halo_tot = np.array(f['a'])
		Rh = halo_tot[3,:]
		rho_g = halo_tot[2,:]
		rho_s = halo_tot[1,:]
		rho_d = halo_tot[0,:]

		rh.append(Rh)
		rhog.append(rho_g)
		rhos.append(rho_s)
		rhod.append(rho_d)

		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/halo_bcg_icm_profile_%s_z_%.3f.h5' % (sub_set[k], z0)) as f:
			cord_array = np.array(f['a'])
		rho_ccm = cord_array[0,:]
		R_ccm = cord_array[1,:]

		rccm.append(R_ccm)
		rhoccm.append(rho_ccm)

		eta_ccm = rho_ccm[rho_d != 0] / rho_d[rho_d != 0]
		etaccm.append(eta_ccm)

	plt.figure(figsize = (16, 12))
	gs = gridspec.GridSpec(2, 1, height_ratios = [3, 2])
	ax = plt.subplot(gs[0])
	bx = plt.subplot(gs[1])
	for kk in range(2):
		ax.plot(rh[kk], rhod[kk], ls = '-', color = mpl.cm.rainbow(kk/2), 
			label = r'$\rho_{DM} \, cluster \, %s$' % sub_set[kk])
		ax.plot(rccm[kk], rhoccm[kk], ls = '--', color = mpl.cm.rainbow(kk/2), 
			label = r'$\rho_{ICM+BCG} \, cluster \, %s$' % sub_set[kk])

		bx.plot(rccm[kk][rhod[kk] != 0], etaccm[kk], ls = '-', color = mpl.cm.rainbow(kk/2), 
			label = r'$\eta_{\rho_{ICM+BCG} / \rho_{DM} \, %s}$' % sub_set[kk])

	ax.set_xlim(1e-1, 2*np.max(rh[1]))
	ax.set_xscale('log')
	ax.set_xlabel(r'$R[kpc/h]$')
	ax.set_ylabel(r'$\rho[h^2 M_\odot /kpc^3]$')
	ax.set_yscale('log')
	ax.legend(loc = 1, fontsize = 10)

	bx.set_xscale('log')
	bx.set_yscale('log')
	bx.set_xlabel(r'$R[arcsec]$')
	bx.set_ylabel(r'$\eta_{\rho}$')
	plt.sca(bx)
	plt.legend(loc = 1, fontsize = 15)

	plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/halo_dens_z0_%s_%s.png' % (sub_set[0], sub_set[1]))
	plt.close()
	raise
	return

def main():
	rho_GX()
	#fig_out()
if __name__ == '__main__':
	main()