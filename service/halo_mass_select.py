import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob
import h5py
import numpy as np
import pandas as pd
import astropy.io.ascii as asc
import astropy.units as U
import astropy.constants as C
from astropy.cosmology import Planck15 as Plank
from scipy.interpolate import interp1d as interp

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

def ICM_GX_subr(exr):
	enl_r = exr * 1
	load2 = 'NewMDCLUSTER_'
	sub_set = CID
	#sub_set = ['0135']
	Nh = len(sub_set)

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

		m_s = np.zeros(N, dtype = np.float)
		m_d = np.zeros(N, dtype = np.float)
		z_n = np.zeros(N, dtype = np.float)

		m_sat = np.zeros(N, dtype = np.float)
		m_bcg = np.zeros(N, dtype = np.float)
		m_ICL = np.zeros(N, dtype = np.float)
		m_ICM = np.zeros(N, dtype = np.float)

		m_s_iner = np.zeros(N, dtype = np.float)
		m_s_out = np.zeros(N, dtype = np.float)
		for q in range(N):
			subz = zc[q]
			subid = snap[N-q]
			if subz <= id_z :
				snap_N = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], subid), 'npartTotal')
				try:
					snap_shot_star = pygdr.readsnap(
					'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], subid), 'pos', 'star')
					snap_mass_star = pygdr.readsnap(
					'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], subid), 'mass','star')
				except SystemExit:
					print('no star particles now')
					break

				halo = asc.read(
				'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/GadgetX-NewMDCLUSTER_%s.snap_%s.z%.3f.AHF_halos'
				 % (sub_set[k], sub_set[k], subid, subz),converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
				Mhalo = np.array(halo['col4'])
				Mgas = np.array(halo['col45'])
				Mstar= np.array(halo['col65'])
				Rvir = np.array(halo['col12'])
				xhalo = np.array(halo['col6'])
				yhalo = np.array(halo['col7'])
				zhalo = np.array(halo['col8'])

				halo_list = np.array(halo['col1'])
				sub_halo = np.array(halo['col2'])

				goal_halo = HID[trc_id, q]
				check_id = find.find1d(halo_list, goal_halo)
				try:
					m_d[q] = Mhalo[check_id]
					m_s[q] = Mstar[check_id]
					z_n[q] = subz * 1
					x0 = xhalo[check_id]
					y0 = yhalo[check_id]
					z0 = zhalo[check_id]
					R0 = Rvir[check_id]

					if snap_N[-2] != 0:
						dstar = np.sqrt((snap_shot_star[:, 0] - x0)**2 + 
							(snap_shot_star[:, 1] - y0)**2 + (snap_shot_star[:, 2] - z0)**2)
						ids = dstar <= R0
						inlstar = snap_shot_star[ids, :]
						inlmass_star = snap_mass_star[ids]
					m_s[q] = np.sum(inlmass_star) * 10**10

					R_range = 150.0
					edge_x = np.array([x0-R_range, x0+R_range])
					edge_y = np.array([y0-R_range, y0+R_range])
					edge_z = np.array([z0-R_range, z0+R_range])
					iv_s = ((inlstar[:,0] <= edge_x[1])&(inlstar[:,0] >= edge_x[0]))&(
							(inlstar[:,1] <= edge_y[1])&(inlstar[:,1] >= edge_y[0]))&(
							(inlstar[:,2] <= edge_z[1])&(inlstar[:,2] >= edge_z[0]))
					inl_star = inlstar[iv_s,:]
					num_bins = np.ceil(R_range*2/resolution)
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
							cenxstar = (bin_x_star[is_max[0]+1] + bin_x_star[is_max[0]])/2.0
							cenystar = (bin_y_star[is_max[1]+1] + bin_y_star[is_max[1]])/2.0
							cenzstar = (bin_z_star[is_max[2]+1] + bin_z_star[is_max[2]])/2.0
							cen_po_star[p,:] = np.array([cenxstar,cenystar,cenzstar])
							hist_use[is_max] = 0.0
						compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+
							(cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
						ismin = find.find1d(compare_d,np.min(compare_d))
						cen_x_star = cen_po_star[ismin,0]
						cen_y_star = cen_po_star[ismin,1]
						cen_z_star = cen_po_star[ismin,2]
						# calculate the bcg mass
						ds_star = np.sqrt((inlstar[:,0] - cen_x_star)**2 + 
							(inlstar[:,1] - cen_y_star)**2 + (inlstar[:,2] - cen_z_star)**2)
						icen = ds_star <= BCG_size
						gx = np.where(icen == True)[0]
						m_bcg[q] = np.sum(inlmass_star[gx]) * 10**10

						# select the satellite and ICL
						try:
							ih = sub_halo == halo_list[check_id]
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
							cenr = real_sub_r[ucen[0]]
							dcr = np.sqrt((inlstar[:,0] - real_sub_x[ucen[0]])**2 + 
								(inlstar[:,1] - real_sub_y[ucen[0]])**2 + (inlstar[:,2] - real_sub_z[ucen[0]])**2)
							cen_mass = np.sum(inlmass_star[dcr <= enl_r * cenr]) * 10**10

							ddm = 0
							for tt in range(len(sub_r)):
								dr = np.sqrt((inlstar[:,0] - sub_x[tt])**2 + 
									(inlstar[:,1] - sub_y[tt])**2 + (inlstar[:,2] - sub_z[tt])**2)
								iic = dr <= enl_r * sub_r[tt]
								ddm = ddm + np.sum(inlmass_star[iic]) * 10**10

							m_sat[q] = ddm - cen_mass
							m_ICL[q] = m_s[q] - m_sat[q] - m_bcg[q]
						except ValueError:
							print('there is noly one halo!')
							print('Now redshift is %.3f' % subz)
							m_sat[q] = 0
							m_ICL[q] = m_s[q] - m_sat[q] - m_bcg[q]

						# for star
						dr_sh = np.sqrt((sub_x - cen_x_star)**2 + (sub_y - cen_y_star)**2 +(sub_z - cen_z_star)**2)
						id_cen = np.where(dr_sh == np.min(dr_sh))[0]

						drs = np.sqrt((inlstar[:, 0] - cen_x_star)**2 +
							(inlstar[:, 1] - cen_y_star)**2 +(inlstar[:, 2] - cen_z_star)**2)
						isr1 = drs <= R_in
						inlstar1 = inlstar[isr1, :]
						inlmass_star1 = inlmass_star[isr1]
						m_s_iner[q] = np.sum(inlmass_star[isr1])*10**10

						isr2 = (drs > R_in) & (drs <= R0)
						inlstar2 = inlstar[isr2, :]
						inlmass_star2 = inlmass_star[isr2]
						m_s_out[q] = np.sum(inlmass_star[isr2])*10**10

						pos = inlstar2 * 1
						out_star = inlmass_star2 * 1

						for tt in range(len(sub_r)):
							if tt == id_cen[0]:
								continue
							else:
								dr = np.sqrt((pos[:,0] - sub_x[tt])**2 + (pos[:,1] - sub_y[tt])**2 + (pos[:,2] - sub_z[tt])**2)
								iic = dr >= enl_r * sub_r[tt]
								pos = pos[iic, :]
								out_star = out_star[iic]
						m_ICM[q] = np.sum(out_star)*10**10
					except ValueError:
						continue
				except IndexError:
					continue
			else:
				continue

		il = m_s != 0
		zgx = np.array(z_n[il])
		ms = np.array(m_s[il])
		md = np.array(m_d[il])

		ms_iner = np.array(m_s_iner[il])
		ms_out = np.array(m_s_out[il])

		ms_sat = np.array(m_sat[il])
		ms_ICL = np.array(m_ICL[il])
		ms_bcg = np.array(m_bcg[il])
		ms_ICM = np.array(m_ICM[il])
		region_data = np.array([zgx, ms, md, ms_iner, ms_out, ms_sat, ms_ICL, ms_bcg, ms_ICM])
		with h5py.File(
				'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/C_%s_mass_trck_with_ICM_%.1fsubr.h5' % (sub_set[k], enl_r), 'w') as f:
			f['a'] = np.array(region_data)
		with h5py.File(
				'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/C_%s_mass_trck_with_ICM_%.1fsubr.h5' % (sub_set[k], enl_r)) as f:
			for kk in range(len(region_data)):
				f['a'][kk, :] = region_data[kk, :]

		print(k)
	return

def rho_GX_subr(exr):
	enl_r = exr * 1
	load2 = 'NewMDCLUSTER_'
	sub_set = CID
	#sub_set = ['0135']

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

		rho_bcg = np.zeros(bins, dtype = np.float)
		R_bcg = np.zeros(bins, dtype = np.float)
		# BCG + ICM component
		rho_ccm = np.zeros((N, 3*bins), dtype = np.float)
		R_ccm = np.zeros((N, 3*bins), dtype = np.float)

		for qq in range(N):
			zq = qq
			z_0 = zc[zq]
			Goal_ID = snap[N - zq]
			if z_0 <= 1:
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
						r_bcg = np.logspace(-1, np.log10(R_BCG), bins + 1)
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
							ih = host_ID == halo_colum[check_id]
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

						except ValueError:
							print('there is noly one halo!')
						# select ICL + BCG from all star 
						drs = np.sqrt((sub_x - cen_x_star)**2 + (sub_y - cen_y_star)**2 +(sub_z - cen_z_star)**2)
						id_cen = np.where(drs == np.min(drs))[0]

						pos = inlstar * 1
						out_star = inlmass_star * 1
						r_icm = np.logspace(-1, np.log10(R0), 3*bins + 1)
						for tt in range(len(sub_r)):
							if tt == id_cen[0]:
								continue
							else:
								dr = np.sqrt((pos[:,0] - sub_x[tt])**2 + (pos[:,1] - sub_y[tt])**2 + (pos[:,2] - sub_z[tt])**2)
								iic = dr >= enl_r * sub_r[tt]
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
							n_subV = 0
							cc_dr = np.sqrt((sub_x - x0)**2 + (sub_y - y0)**2 + (sub_z - z0)**2)
							for ss in range(len(sub_r)):
								if tt == id_cen[0]:
									V_sub = 0
									continue
								elif (((r_icm[kk+1] >= cc_dr[ss] - enl_r * sub_r[ss]) & (r_icm[kk+1] <= cc_dr[ss] + enl_r * sub_r[ss])) | 
									((r_icm[kk] >= cc_dr[ss] - enl_r * sub_r[ss]) & (r_icm[kk] <= enl_r * cc_dr[ss] + sub_r[ss]))):
									cen_ds0 = np.abs(cc_dr[ss] - r_icm[kk])
									cen_ds1 = np.abs(cc_dr[ss] - r_icm[kk+1])
									if cen_ds0 > enl_r * sub_r[ss]:
										if np.abs(cen_ds0 - enl_r * sub_r[ss]) >= 0.5 * ddr:
											s1 = np.pi * (np.sqrt(enl_r * sub_r[ss]**2 - cen_ds1**2))**2
											s0 = 0.5 * s1
											V_sub = (s1 + s0 + np.sqrt(s0 * s1)) * ddr / 3
										else:
											V_sub = 0
									elif cen_ds1 > enl_r * sub_r[ss]:
										if np.abs(enl_r * sub_r[ss] - cen_ds0) >= 0.5 * ddr:
											s0 = np.pi * (np.sqrt(enl_r * sub_r[ss]**2 - cen_ds0**2))**2
											s1 = 0.5 * s0
											V_sub = (s0 + s1 + np.sqrt(s0 * s1)) * ddr / 3
										else:
											V_sub = 0
									else:
										s0 = np.pi * (np.sqrt(enl_r * sub_r[ss]**2 - cen_ds0**2))**2
										s1 = np.pi * (np.sqrt(enl_r * sub_r[ss]**2 - cen_ds1**2))**2
										V_sub = (s0 + s1 + np.sqrt(s0 * s1)) * ddr / 3
								else:
									V_sub = 0
									continue
								n_subV = n_subV + V_sub
							rho_ccm[qq, kk] = (dm_icm[kk + 1] - dm_icm[kk]) / (4*np.pi * ddr * r_icm[kk]**2 - n_subV)
							R_ccm[qq, kk] = 0.5 * (r_icm[kk + 1] + r_icm[kk])
					except ValueError:
						continue
				except IndexError:
					continue
		cord_array = np.array([rho_ccm, R_ccm])
		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/Cluster_profile/halo_bcg_icm_profile_%s_%.1fsubr.h5' % (sub_set[k], enl_r), 'w') as f:
			f['a'] = np.array(cord_array)
		with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/Cluster_profile/halo_bcg_icm_profile_%s_%.1fsubr.h5' % (sub_set[k], enl_r)) as f:
			for qq in range(len(cord_array)):
				f['a'][qq,:] = cord_array[qq,:]
		print(k)
	return

def fig_out(exr):
	exr = exr
	eps = 1e-2
	# mass accrection
	#sub_set = ['0135']
	sub_set = CID
	Nh = len(sub_set)

	for k in range(Nh):
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/C_%s_mass_trck_with_ICM.h5' % sub_set[k]) as f:
			record = np.array(f['a'])
		z0 = record[0,:]
		Ms0 = record[1,:]
		Mh0 = record[2,:]
		Msat0 = record[5,:]
		Mbcg0 = record[7,:]
		Micm0 = record[8,:]
		t_Ms0 = Msat0 + Mbcg0 + Micm0

		eta_sat0 = Msat0[Msat0 >= 10**11] / Mh0[Msat0 >= 10**11]
		zsat0 = z0[Msat0 >= 10**11]
		eta_bcg0 = Mbcg0[Mbcg0 >= 10**11] / Mh0[Mbcg0 >= 10**11]
		zbcg0 = z0[Mbcg0 >= 10**11]
		eta_icm0 = Micm0[Micm0 >= 10**11] / Mh0[Micm0 >= 10**11]
		zicm0 = z0[Micm0 >= 10**11]
		eta_tot_m0 = t_Ms0[t_Ms0 >= 10**11] / Mh0[t_Ms0 >= 10**11]
		z_tot0 = z0[t_Ms0 >= 10**11]

		sub_z = []
		sub_ms = []
		sub_mh = []
		sub_sat = []
		sub_bcg = []
		sub_icm = []
		sub_tot = []

		sub_e_sat = []
		sub_e_zsat = []
		sub_e_bcg = []
		sub_e_zbcg = []
		sub_e_icm = []
		sub_e_zicm = []
		sub_e_tot = []
		sub_e_ztot = []
		for t in range(len(exr)):
			enl_r = exr[t]
			with h5py.File(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/C_%s_mass_trck_with_ICM_%.1fsubr.h5' % (sub_set[k], enl_r)) as f:
				sub_data = np.array(f['a'])
			zgx = sub_data[0,:]
			ms = sub_data[1,:]
			mh = sub_data[2,:]
			msat = sub_data[5,:]
			mbcg = sub_data[7,:]
			micm = sub_data[8,:]
			mtot = msat + mbcg + micm

			eta_sat = msat[msat >= 10**11] / mh[msat >= 10**11]
			zsat = zgx[msat >= 10**11]
			eta_bcg = mbcg[mbcg >= 10**11] / mh[mbcg >= 10**11]
			zbcg = zgx[mbcg >= 10**11]
			eta_icm = micm[micm >= 10**11] / mh[micm >= 10**11]
			zicm = zgx[micm >= 10**11]
			eta_tot_m = mtot[mtot >= 10**11] / mh[mtot >= 10**11]
			z_tot = zgx[mtot >= 10**11]

			sub_z.append(zgx)
			sub_ms.append(ms)
			sub_mh.append(mh)
			sub_sat.append(msat)
			sub_bcg.append(mbcg)
			sub_icm.append(micm)
			sub_tot.append(mtot)

			sub_e_sat.append(eta_sat)
			sub_e_zsat.append(zsat)
			sub_e_bcg.append(eta_bcg)
			sub_e_zbcg.append(zbcg)
			sub_e_icm.append(eta_icm)
			sub_e_zicm.append(zicm)
			sub_e_tot.append(eta_tot_m)
			sub_e_ztot.append(z_tot)
		## fig part
		plt.figure(figsize = (12, 12))
		gs = gridspec.GridSpec(2, 1, height_ratios = [3, 2])
		ax = plt.subplot(gs[0])
		bx = plt.subplot(gs[1], sharex = ax)

		ax.set_title('Assembly history', fontsize = 15)
		ax.plot(np.log10(1 + z0[(Mh0 != 0) & (Mh0* eps >= 10**11)]),
			Mh0[(Mh0 != 0) & (Mh0* eps >= 10**11)] * eps, 'k-', label = r'$M_{h}/10^{2} \, R_{vir}$')
		ax.plot(np.log10(1 + z0[(Msat0 != 0) & (Msat0 >= 10**11)]),
			Msat0[(Msat0 != 0) & (Msat0 >= 10**11)], 'b-', label = r'$M_{sat} \, subR_{vir}$')
		bx.plot(np.log10(1 + zsat0), eta_sat0, 'b-', label = r'$M_{sat} / M_h \, subR_{vir}$')
		for tt in range(len(exr)):
			ax.plot(np.log10(1 + sub_z[tt][sub_sat[tt] >= 10**11]), sub_sat[tt][sub_sat[tt] >= 10**11], 
				ls = '-', color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{sat} \, %.1f subR_{vir}$' % exr[tt])
			bx.plot(np.log10(1 + sub_e_zsat[tt]), sub_e_sat[tt], ls = '-', 
				color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{sat} / M_h \, %.1f subR_{vir}$' % exr[tt])
			# handles, labels = plt.gca().get_legend_handles_labels()
		ax1 = ax.twiny()
		xtik = ax.get_xticks()
		Zt = 10**(xtik) - 1
		LBT = Plank.lookback_time(Zt).value
		ax1.set_xticks(xtik)
		ax1.set_xticklabels(["%.2f" % ll for ll in LBT])
		ax1.set_xlim(ax.get_xlim())
		ax1.set_xlabel(r'$look \; back \; time[Gyr]$', fontsize = 12)

		ax.set_ylabel(r'$M[M_{\odot}/h]$', fontsize = 12)
		ax.set_yscale('log')
		ax.legend(loc = 1, fontsize = 12)
		ax.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		ax1.tick_params(axis = 'x', which = 'both', direction = 'in', labelsize = 12)

		bx.set_xlabel(r'$log(1+z)$', fontsize = 12)
		bx.set_ylabel(r'$Mass \; ratio$', fontsize = 12)
		bx.set_yscale('log')
		bx.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		bx.legend(loc = 4, fontsize = 12)

		plt.subplots_adjust(hspace = 0)
		plt.savefig(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/compare_fig/C_%s_Msat_accret_with_differ_subr.png' % sub_set[k], dpi = 300)
		plt.close()

		plt.figure(figsize = (12, 12))
		gs = gridspec.GridSpec(2, 1, height_ratios = [3, 2])
		ax = plt.subplot(gs[0])
		bx = plt.subplot(gs[1], sharex = ax)
		ax.plot(np.log10(1 + z0[(Mh0 != 0) & (Mh0* eps >= 10**11)]),
			Mh0[(Mh0 != 0) & (Mh0* eps >= 10**11)] * eps, 'k-', label = r'$M_{h}/10^{2} \, subR_{vir}$')
		ax.plot(np.log10(1 + z0[(Mbcg0 != 0) & (Mbcg0 >= 10**11)]),
			Mbcg0[(Mbcg0 != 0) & (Mbcg0 >= 10**11)], 'b-', label = r'$M_{BCG} \, subR_{vir}$')
		bx.plot(np.log10(1 + zbcg0), eta_bcg0, 'b-', label = r'$M_{BCG} / M_h \, subR_{vir}$')
		for tt in range(len(exr)):
			ax.plot(np.log10(1 + sub_z[tt][sub_bcg[tt] >= 10**11]), sub_bcg[tt][sub_bcg[tt] >= 10**11], 
				ls = '-', color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{BCG} \, %.1f subR_{vir}$' % exr[tt])
			bx.plot(np.log10(1 + sub_e_zbcg[tt]), sub_e_bcg[tt], ls = '-', 
				color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{BCG} / M_h \, %.1f subR_{vir}$' % exr[tt])
			# handles, labels = plt.gca().get_legend_handles_labels()
		ax1 = ax.twiny()
		xtik = ax.get_xticks()
		Zt = 10**(xtik) - 1
		LBT = Plank.lookback_time(Zt).value
		ax1.set_xticks(xtik)
		ax1.set_xticklabels(["%.2f" % ll for ll in LBT])
		ax1.set_xlim(ax.get_xlim())
		ax1.set_xlabel(r'$look \; back \; time[Gyr]$', fontsize = 12)

		ax.set_ylabel(r'$M[M_{\odot}/h]$', fontsize = 12)
		ax.set_yscale('log')
		ax.legend(loc = 1, fontsize = 12)
		ax.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		ax1.tick_params(axis = 'x', which = 'both', direction = 'in', labelsize = 12)

		bx.set_xlabel(r'$log(1+z)$', fontsize = 12)
		bx.set_ylabel(r'$Mass \; ratio$', fontsize = 12)
		bx.set_yscale('log')
		bx.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		bx.legend(loc = 4, fontsize = 12)

		plt.subplots_adjust(hspace = 0)
		plt.savefig(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/compare_fig/C_%s_Mbcg_accret_with_differ_subr.png' % sub_set[k], dpi = 300)
		plt.close()

		plt.figure(figsize = (12, 12))
		gs = gridspec.GridSpec(2, 1, height_ratios = [3, 2])
		ax = plt.subplot(gs[0])
		bx = plt.subplot(gs[1], sharex = ax)
		ax.plot(np.log10(1 + z0[(Mh0 != 0) & (Mh0* eps >= 10**11)]),
			Mh0[(Mh0 != 0) & (Mh0* eps >= 10**11)] * eps, 'k-', label = r'$M_{h}/10^{2} \, subR_{vir}$')
		ax.plot(np.log10(1 + z0[(Micm0 != 0) & (Micm0 >= 10**11)]), 
			Micm0[(Micm0 != 0) & (Micm0 >= 10**11)], 'b-', label = r'$M_{ICM} \, subR_{vir}$')
		bx.plot(np.log10(1 + zicm0), eta_icm0, 'b-', label = r'$M_{ICM} / M_h \, subR_{vir}$')
		for tt in range(len(exr)):
			ax.plot(np.log10(1 + sub_z[tt][sub_icm[tt] >= 10**11]), sub_icm[tt][sub_icm[tt] >= 10**11], 
				ls = '-', color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{ICM} \, %.1f subR_{vir}$' % exr[tt])
			bx.plot(np.log10(1 + sub_e_zicm[tt]), sub_e_icm[tt], ls = '-', 
				color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{ICM} / M_h \, %.1f subR_{vir}$' % exr[tt])
			# handles, labels = plt.gca().get_legend_handles_labels()
		ax1 = ax.twiny()
		xtik = ax.get_xticks()
		Zt = 10**(xtik) - 1
		LBT = Plank.lookback_time(Zt).value
		ax1.set_xticks(xtik)
		ax1.set_xticklabels(["%.2f" % ll for ll in LBT])
		ax1.set_xlim(ax.get_xlim())
		ax1.set_xlabel(r'$look \; back \; time[Gyr]$', fontsize = 12)

		ax.set_ylabel(r'$M[M_{\odot}/h]$', fontsize = 12)
		ax.set_yscale('log')
		ax.legend(loc = 1, fontsize = 12)
		ax.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		ax1.tick_params(axis = 'x', which = 'both', direction = 'in', labelsize = 12)

		bx.set_xlabel(r'$log(1+z)$', fontsize = 12)
		bx.set_ylabel(r'$Mass \; ratio$', fontsize = 12)
		bx.set_yscale('log')
		bx.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		bx.legend(loc = 4, fontsize = 12)

		plt.subplots_adjust(hspace = 0)
		plt.savefig(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/compare_fig/C_%s_Micm_accret_with_differ_subr.png' % sub_set[k], dpi = 300)
		plt.close()

		plt.figure(figsize = (12, 12))
		gs = gridspec.GridSpec(2, 1, height_ratios = [3, 2])
		ax = plt.subplot(gs[0])
		bx = plt.subplot(gs[1], sharex = ax)
		ax.plot(np.log10(1 + z0[(Mh0 != 0) & (Mh0* eps >= 10**11)]),
			Mh0[(Mh0 != 0) & (Mh0* eps >= 10**11)] * eps, 'k-', label = r'$M_{h}/10^{2} \, subR_{vir}$')		
		ax.plot(np.log10(1 + z0[t_Ms0 >= 10**11]), t_Ms0[t_Ms0 >= 10**11], 'b-', label = r'$M^{\ast}_{BCG + ICM + sat} \, subR_{vir}$')
		bx.plot(np.log10(1 + z_tot0), eta_tot_m0, 'b-', label = r'$M_{BCG + ICM + sat} / M_h \, subR_{vir}$')
		for tt in range(len(exr)):
			ax.plot(np.log10(1 + sub_z[tt][sub_tot[tt] >= 10**11]), sub_tot[tt][sub_tot[tt] >= 10**11], 
				ls = '-', color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{BCG + ICM + sat} \, %.1f subR_{vir}$' % exr[tt])
			bx.plot(np.log10(1 + sub_e_ztot[tt]), sub_e_tot[tt], ls = '-', 
				color = mpl.cm.plasma(tt/len(exr)), label = r'$M_{BCG + ICM + sat} / M_h \, %.1f subR_{vir}$' % exr[tt])
			# handles, labels = plt.gca().get_legend_handles_labels()
		ax1 = ax.twiny()
		xtik = ax.get_xticks()
		Zt = 10**(xtik) - 1
		LBT = Plank.lookback_time(Zt).value
		ax1.set_xticks(xtik)
		ax1.set_xticklabels(["%.2f" % ll for ll in LBT])
		ax1.set_xlim(ax.get_xlim())
		ax1.set_xlabel(r'$look \; back \; time[Gyr]$', fontsize = 12)

		ax.set_ylabel(r'$M[M_{\odot}/h]$', fontsize = 12)
		ax.set_yscale('log')
		ax.legend(loc = 1, fontsize = 12)
		ax.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		ax1.tick_params(axis = 'x', which = 'both', direction = 'in', labelsize = 12)

		bx.set_xlabel(r'$log(1+z)$', fontsize = 12)
		bx.set_ylabel(r'$Mass \; ratio$', fontsize = 12)
		bx.set_yscale('log')
		bx.tick_params(axis = 'both', which = 'both', direction = 'in', labelsize = 12)
		bx.legend(loc = 4, fontsize = 12)

		plt.subplots_adjust(hspace = 0)
		plt.savefig(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/compare_fig/C_%s_Mtot_accret_with_differ_subr.png' % sub_set[k], dpi = 300)
		plt.close()
	raise
	return

def main():
	ex_r = np.linspace(1.2, 1.6, 3)
	for pp in range(len(ex_r)):
		ICM_GX_subr(ex_r[pp])
		rho_GX_subr(ex_r[pp])
	fig_out(ex_r)

if __name__ == '__main__':
	main()
