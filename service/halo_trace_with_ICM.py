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

def ICM_GX():

	load2 = 'NewMDCLUSTER_'
	Nh = len(CID)

	for k in range(Nh):

		if (CID[k] == '0080') | (CID[k] == '0110'):
			trc_id = 1
		else:
			trc_id = 0

		sub_file = glob.glob(
			'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/*.AHF_halos' % CID[k])
		snap_id = [str(sub_file[107:110]) for sub_file in sub_file]
		snap = np.sort(snap_id)

		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/z_for_C_%s.h5' % CID[k]) as f:
			zc = np.array(f['a'])

		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/cluster_%s_tree.h5' % CID[k]) as f:
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
		m_s_media = np.zeros(N, dtype = np.float)
		m_s_out = np.zeros(N, dtype = np.float)
		for q in range(N):
			subz = zc[q]
			subid = snap[N-q]
			if subz <= id_z :
				snap_N = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (CID[k], subid), 'npartTotal')
				try:
					snap_shot_star = pygdr.readsnap(
					'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (CID[k], subid), 'pos', 'star')
					snap_mass_star = pygdr.readsnap(
					'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (CID[k], subid), 'mass','star')
				except SystemExit:
					print('no star particles now')
					break

				halo = asc.read(
				'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/GadgetX-NewMDCLUSTER_%s.snap_%s.z%.3f.AHF_halos'
				 % (CID[k], CID[k], subid, subz),converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
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
							cen_mass = real_sub[ucen[0]]
							m_sat[q] = np.sum(real_sub) - cen_mass
							m_ICL[q] = m_s[q] - m_sat[q] - m_bcg[q]
						except ValueError:
							print('there is noly one halo!')
							print('Now redshift is %.3f' % subz)
							m_sat[q] = 0
							m_ICL[q] = m_s[q] - m_sat[q] - m_bcg[q]

						# for star 
						drs = np.sqrt((inlstar[:, 0] - cen_x_star)**2 +
							(inlstar[:, 1] - cen_y_star)**2 +(inlstar[:, 2] - cen_z_star)**2)
						isr1 = drs <= R_in
						inlstar1 = inlstar[isr1, :]
						inlmass_star1 = inlmass_star[isr1]
						m_s_iner[q] = np.sum(inlmass_star[isr1])*10**10

						isr2 = (drs > R_in) & (drs <= R_out)
						inlstar2 = inlstar[isr2, :]
						inlmass_star2 = inlmass_star[isr2]
						m_s_media[q] = np.sum(inlmass_star[isr2])*10**10

						isr3 = (drs > R_out) & (drs <= R0)
						inlstar3 = inlstar[isr3, :]
						inlmass_star3 = inlmass_star[isr3]
						m_s_out[q] = np.sum(inlmass_star[isr3])*10**10

						pos = inlstar3*1
						out_star = inlmass_star3*1
						for tt in range(len(sub_r)):
							dr = np.sqrt((pos[:,0] - sub_x[tt])**2 + 
								(pos[:,1] - sub_y[tt])**2 + (pos[:,2] - sub_z[tt])**2)
							iic = dr >= sub_r[tt]
							pos = pos[iic, :]
							out_star = out_star[iic]
						m_ICM[q] = np.sum(out_star)*10**10
					except ValueError:
						print('There is no enough star to devide bins to find BCG,redshift is %.3f'%snap_z)
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
		ms_media = np.array(m_s_media[il])
		ms_out = np.array(m_s_out[il])

		ms_sat = np.array(m_sat[il])
		ms_ICL = np.array(m_ICL[il])
		ms_bcg = np.array(m_bcg[il])
		ms_ICM = np.array(m_ICM[il])
		region_data = np.array([zgx, ms, md, 
						ms_iner, ms_media, ms_out, 
						ms_sat, ms_ICL, ms_bcg, ms_ICM])
		with h5py.File(
				'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/C_%s_mass_trck_with_ICM.h5' % CID[k], 'w') as f:
			f['a'] = np.array(region_data)
		with h5py.File(
				'/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/C_%s_mass_trck_with_ICM.h5' % CID[k]) as f:
			for kk in range(len(region_data)):
				f['a'][kk, :] = region_data[kk, :]

		print(k)
	return

def fig_out():
	#sub_set = ['0145', '0182'] 
	sub_set = ['0153', '0182']
	setN = len(sub_set)
	Z = []
	M_sat = []
	M_bcg = []
	M_icm = []
	M_h = []
	for k in range(setN):
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/C_%s_mass_trck_with_ICM.h5' % sub_set[k]) as f:
			record = np.array(f['a'])
		zgx = record[0,:]
		mh = record[2,:]
		msat = record[6,:]
		mbcg = record[8,:]
		micm = record[9,:]

		Z.append(zgx)
		M_h.append(mh)
		M_sat.append(msat)
		M_bcg.append(mbcg)
		M_icm.append(micm)

	Zt = Z[1]
	tL = Plank.lookback_time(Zt)

	zh0 = np.array(Z[0])
	zh1 = np.array(Z[1])
	Mh0 = np.array(M_h[0])
	Mh1 = np.array(M_h[1])
	eta_dm0 = Mh0 / Mh0[0]
	eta_dm1 = Mh1 / Mh1[0]

	Msat0 = np.array(M_sat[0])
	Msat1 = np.array(M_sat[1])
	
	eta_sat0 = Msat0[Msat0 != 0] / Mh0[Msat0 != 0]
	eta_sat1 = Msat1[Msat1 != 0] / Mh1[Msat1 != 0]
	'''
	eta_sat0 = Msat0[Msat0 != 0] / Mh0[0]
	eta_sat1 = Msat1[Msat1 != 0] / Mh1[0]
	'''
	
	Mbcg0 = np.array(M_bcg[0])
	Mbcg1 = np.array(M_bcg[1])
	
	eta_bcg0 = Mbcg0[Mbcg0 != 0] / Mh0[Mbcg0 != 0]
	eta_bcg1 = Mbcg1[Mbcg1 != 0] / Mh1[Mbcg1 != 0]
	'''
	eta_bcg0 = Mbcg0[Mbcg0 != 0] / Mh0[0]
	eta_bcg1 = Mbcg1[Mbcg1 != 0] / Mh1[0]
	'''

	Micm0 = np.array(M_icm[0])
	Micm1 = np.array(M_icm[1])
	
	eta_icm0 = Micm0[Micm0 != 0] / Mh0[Micm0 != 0]
	eta_icm1 = Micm1[Micm1 != 0] / Mh1[Micm1 != 0]
	'''
	eta_icm0 = Micm0[Micm0 != 0] / Mh0[0]
	eta_icm1 = Micm1[Micm1 != 0] / Mh1[0]
	'''

	plt.figure(figsize = (16, 12))
	gs = gridspec.GridSpec(2, 1, height_ratios = [3, 2])
	ax = plt.subplot(gs[0])
	bx = plt.subplot(gs[1])

	ax.set_title('Assembly comparation')
	ax.plot(np.log10(1 + zh0[Mh0 != 0]), Mh0[Mh0 != 0], 'k-', label = r'$M_{h} \, C_{%s}$' % sub_set[0])
	ax.plot(np.log10(1 + zh0[Msat0 != 0]), Msat0[Msat0 != 0], 'b-', label = r'$M_{sat} \, C_{%s}$' % sub_set[0])
	ax.plot(np.log10(1 + zh0[Mbcg0 != 0]), Mbcg0[Mbcg0 != 0], 'r-', label = r'$M_{BCG} \, C_{%s}$' % sub_set[0])
	ax.plot(np.log10(1 + zh0[Micm0 != 0]), Micm0[Micm0 != 0], 'g-', label = r'$M_{ICM} \, C_{%s}$' % sub_set[0])
	ax.plot(np.log10(1 + zh1[Mh1 != 0]), Mh1[Mh1 != 0], 'k--', label = r'$M_{h} \, C_{%s}$' % sub_set[1])
	ax.plot(np.log10(1 + zh1[Msat1 != 0]), Msat1[Msat1 != 0], 'b--', label = r'$M_{sat} \, C_{%s}$' % sub_set[1])
	ax.plot(np.log10(1 + zh1[Mbcg1 != 0]), Mbcg1[Mbcg1 != 0], 'r--', label = r'$M_{BCG} \, C_{%s}$' % sub_set[1])
	ax.plot(np.log10(1 + zh1[Micm1 != 0]), Micm1[Micm1 != 0], 'g--', label = r'$M_{ICM} \, C_{%s}$' % sub_set[1])

	ax1 = ax.twiny()
	ax1.plot(tL[Mh1 != 0], Mh1[Mh1 != 0], 'w--')
	#ax1.set_xscale('log')
	ax1.set_xlabel(r'$look \; back \; time[Gyr]$')
	ax1.set_yscale('log')

	ax.set_xlabel('log(1+z)')
	ax.set_ylabel(r'$M[M_{\odot}/h]$')
	ax.set_yscale('log')
	#ax.set_xscale('log')
	ax.legend(loc = 1)
	ax.tick_params(axis = 'both', which = 'both', direction = 'in')
	ax1.tick_params(axis = 'x', which = 'both', direction = 'in')
	
	bx.plot(np.log10(1 + zh0[Msat0 != 0]), eta_sat0, 'b-', label = r'$\eta_{M_{sat} / M_h} \, %s$' % sub_set[0])
	bx.plot(np.log10(1 + zh0[Mbcg0 != 0]), eta_bcg0, 'r-', label = r'$\eta_{M_{BCG} / M_h} \, %s$' % sub_set[0])
	bx.plot(np.log10(1 + zh0[Micm0 != 0]), eta_icm0, 'g-', label = r'$\eta_{M_{ICM} / M_h} \, %s$' % sub_set[0])
	bx.plot(np.log10(1 + zh1[Msat1 != 0]), eta_sat1, 'b--', label = r'$\eta_{M_{sat} / M_h} \, %s$' % sub_set[1])
	bx.plot(np.log10(1 + zh1[Mbcg1 != 0]), eta_bcg1, 'r--', label = r'$\eta_{M_{BCG} / M_h} \, %s$' % sub_set[1])
	bx.plot(np.log10(1 + zh1[Micm1 != 0]), eta_icm1, 'g--', label = r'$\eta_{M_{ICM} / M_h} \, %s$' % sub_set[1])

	#bx.set_xscale('log')
	bx.set_xlabel(r'$log(1+z)$')
	bx.set_ylabel(r'$\eta_{M}$')
	bx.set_yscale('log')
	bx.tick_params(axis = 'both', which = 'both', direction = 'in')
	bx.legend(loc = 1)

	plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/assembly_comparation_%s_%s.png' % (sub_set[0], sub_set[1]), dpi = 300)
	#plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/assembly_comparation_%s_%s_on_Mh0.png' % (sub_set[0], sub_set[1]), dpi = 300)
	plt.close()
	return

def main():
	#ICM_GX()
	fig_out()

if __name__ == '__main__':
	main()
