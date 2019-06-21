# use for check native mass
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import handy.scatter as hsc

import h5py
import find
import numpy as np
import pygadgetreader as pygdr
import astropy.io.ascii as asc
import scipy.stats as stats
 
resolution = 4.
_id_ = 0 # trace the first one halo(at z = 0)
BCG_size = 50 # in unit kpc/h
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)

with h5py.File('/home/cxkttwl/Scatter/MUSIC/Redshift.h5') as f:
	y0 = f['a']
	com_z = np.array(y0)
with h5py.File('/home/cxkttwl/Scatter/MUSIC/main_tree.h5') as f:
	y1 = f['a']
	tree_line = np.array(y1)
with h5py.File('/home/cxkttwl/Scatter/GadgetX/Redshift_GX.h5') as f:
	y2 = f['a']
	com_z_gx = np.array(y2)
with h5py.File('/home/cxkttwl/Scatter/GadgetX/main_tree_GX.h5') as f:
	y3 = f['a']
	tree_line_gx = np.array(y3)
L = tree_line.shape[0]
iv = np.zeros(L,dtype = np.int0)
for l in range(L):
	u = find.find1d(tree_line[l,:],0)
	iv[l] = u
L_gx = tree_line_gx.shape[0]
iv_gx = np.zeros(L_gx,dtype = np.int0)
for ll in range(L_gx):
	u_gx = find.find1d(tree_line_gx[ll,:],0)
	iv_gx[ll] = u_gx
id_z = com_z[iv[_id_]]
id_z_gx = com_z_gx[iv_gx[_id_]]

for k in range(len(alpha)):
	dd = alpha[-k]
	ss = str(dd)
	if len(ss) == 1:
		s_u = str('00%.0f' % dd)
	elif len(ss) == 2:
		s_u = str('0%.0f' % dd)
	else:
		s_u = str(dd)
	No_snap = s_u
	snap_z = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
	snap_z = np.abs(snap_z)
	print('Now redshift is %.3f'%snap_z)
	snap_N = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
	if ((snap_z >= 2.75) & (snap_z <= 4.3)):
		zt = np.float('%.3f'%snap_z)
		snap_shot_gas = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
		snap_shot_DM = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
		try:
			snap_shot_star = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
		except SystemExit:
			print('no star particles now')
			break
		main_halo = asc.read(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos' % snap_z,
			converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
		goal_z = np.float('%.3f'%snap_z)
		goalz = find.find1d(com_z,goal_z)
		goal_halo = tree_line[_id_,goalz]
		#check the halo in this redshift
		HALO = np.array(main_halo['col1'])
		check_id = find.find1d(HALO,goal_halo) 
		Rvir = np.array(main_halo['col12'])
		xhalo = np.array(main_halo['col6'])
		yhalo = np.array(main_halo['col7'])
		zhalo = np.array(main_halo['col8'])
		x0 = xhalo[check_id]
		y0 = yhalo[check_id]
		z0 = zhalo[check_id]
		R0 = Rvir[check_id]

		dgas = np.sqrt((snap_shot_gas[:, 0]-x0)**2+(snap_shot_gas[:, 1]-y0)**2 +
			(snap_shot_gas[:, 2]-z0)**2)
		ig = dgas <= 1.5*R0
		inl_gas = snap_shot_gas[ig, :]

		dDM = np.sqrt((snap_shot_DM[:, 0]-x0)**2+(snap_shot_DM[:, 1]-y0)**2 +
			(snap_shot_DM[:, 2]-z0)**2)
		iD = dDM <= 1.5*R0
		inl_DM = snap_shot_DM[iD, :]

		if snap_N[-2] != 0:
			dstar = np.sqrt((snap_shot_star[:, 0]-x0)**2+(snap_shot_star[:, 1]-y0)**2 +
				(snap_shot_star[:, 2]-z0)**2)
			ids = dstar <= 1.5*R0
			inl_star = snap_shot_star[ids, :]

		numbins = np.int(np.ceil(2*R0/resolution))
		plt.figure(figsize=(8, 8))
		'''
		plt.hist2d(inl_DM[:, 0], inl_DM[:, 1], bins=[numbins, numbins],
			cmap='viridis', vmin=1e-1, vmax=len(inl_DM)/100, norm=mpl.colors.LogNorm())
		plt.hist2d(inl_gas[:, 0], inl_gas[:, 1], bins=[numbins, numbins],
			cmap='cool', vmin=1e-1, vmax=len(inl_gas)/100, norm=mpl.colors.LogNorm())
		if snap_N[-2] != 0:
			plt.hist2d(inl_star[:, 0], inl_star[:, 1], bins=[numbins, numbins],
				cmap='plasma', vmin = 1e-1, vmax = len(inl_star)/1000, norm=mpl.colors.LogNorm())
		'''
		if snap_N[-2] != 0:
			bind1 = stats.binned_statistic_2d(inl_star[:, 0], inl_star[:, 1], inl_star[:,2],
				statistic = 'count', bins=[numbins, numbins],)
			Nstar = bind1.statistic
			bindx = np.linspace(np.min(inl_star[:, 0]), np.max(inl_star[:, 0]), numbins)
			bindy = np.linspace(np.min(inl_star[:, 1]), np.max(inl_star[:, 1]), numbins)
			plt.pcolormesh(bindx, bindy, Nstar, cmap='rainbow', vmin=1e-1, vmax=np.max(Nstar), 
				alpha=1, norm = mpl.colors.LogNorm())
			plt.plot(x0, y0, '+', c = 'b')
			hsc.circles(x0, y0, s = R0, fc = '', ec = 'r')
			plt.title(r'$GadgetMUSIC 1.5R_{200} z_{%.3f} Re_{%.3f}$' % (snap_z, resolution))
			plt.xlim(x0 - 1.5*R0, x0 + 1.5*R0)
			plt.ylim(y0 - 1.5*R0, y0 + 1.5*R0)
			plt.xlabel('r[kpc/h]')
			plt.ylabel('r[kpc/h]')		
			plt.savefig('MU-z%.3f-distribution.png' % zt, dpi=600)
			plt.show()
			plt.close()

for k in range(len(alpha)):
	dd = alpha[-k]
	ss = str(dd)
	if len(ss) == 1:
		s_u = str('00%.0f' % dd)
	elif len(ss) == 2:
		s_u = str('0%.0f' % dd)
	else:
		s_u = str(dd)
	No_snap = s_u
	snap_z_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
	snap_z_gx = np.abs(snap_z_gx)
	print('Now redshift is %.3f'%snap_z_gx)
	snap_N_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
	if (snap_z_gx >= 2.75) & (snap_z_gx <= 4.3):
		zt = np.float('%.3f'%snap_z_gx)
		snap_shot_gas_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
		snap_shot_DM_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
		try:
			snap_shot_star_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
		except SystemExit:
			print('no star particles now')
			break
		main_halo_gx = asc.read(
			'/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos' % snap_z_gx,
			converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
		goal_z_gx = np.float('%.3f'%snap_z_gx)
		goalz_gx = find.find1d(com_z_gx,goal_z_gx)
		goal_halo_gx = tree_line_gx[_id_,goalz_gx]
		#check the halo in this redshift
		HALO_gx = np.array(main_halo_gx['col1'])
		check_id_gx = find.find1d(HALO_gx,goal_halo_gx) 
		Rvir_gx = np.array(main_halo_gx['col12'])
		xhalo_gx = np.array(main_halo_gx['col6'])
		yhalo_gx = np.array(main_halo_gx['col7'])
		zhalo_gx = np.array(main_halo_gx['col8'])
		x0_gx = xhalo_gx[check_id_gx]
		y0_gx = yhalo_gx[check_id_gx]
		z0_gx = zhalo_gx[check_id_gx]
		R0_gx = Rvir_gx[check_id_gx] 

		dgas_gx = np.sqrt((snap_shot_gas_gx[:, 0]-x0_gx)**2+(snap_shot_gas_gx[:, 1]-y0_gx)**2 +
			(snap_shot_gas_gx[:, 2]-z0_gx)**2)
		ig_gx = dgas_gx <= 1.5*R0_gx
		inl_gas_gx = snap_shot_gas_gx[ig_gx, :]

		dDM_gx = np.sqrt((snap_shot_DM_gx[:, 0]-x0_gx)**2+(snap_shot_DM_gx[:, 1]-y0_gx)**2 +
			(snap_shot_DM_gx[:, 2]-z0_gx)**2)
		iD_gx = dDM_gx <= 1.5*R0_gx
		inl_DM_gx = snap_shot_DM_gx[iD_gx, :]

		if snap_N_gx[-2] != 0:
			dstar_gx = np.sqrt((snap_shot_star_gx[:, 0]-x0_gx)**2+(snap_shot_star_gx[:, 1]-y0_gx)**2 +
				(snap_shot_star_gx[:, 2]-z0_gx)**2)
			ids_gx = dstar_gx <= 1.5*R0_gx
			inl_star_gx = snap_shot_star_gx[ids_gx, :]

		numbins_gx = np.int(np.ceil(2*R0_gx/resolution))
		plt.figure(figsize=(8, 8))
		'''
		plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 1], bins=[numbins_gx, numbins_gx],
			cmap='viridis', vmin=1e-1, vmax=len(inl_DM_gx)/100, norm=mpl.colors.LogNorm())
		plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 1], bins=[numbins_gx, numbins_gx],
			cmap='cool', vmin=1e-1, vmax=len(inl_gas_gx)/100, norm=mpl.colors.LogNorm())
		if snap_N_gx[-2] != 0:
			plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 1], bins=[numbins_gx, numbins_gx],
				cmap='plasma', vmin = 1e-1, vmax = len(inl_star_gx)/1000, norm=mpl.colors.LogNorm())
		'''
		if snap_N_gx[-2] != 0:
			bind2 = stats.binned_statistic_2d(inl_star_gx[:, 0], inl_star_gx[:, 1], inl_star_gx[:,2],
				statistic = 'count', bins=[numbins_gx, numbins_gx],)
			Nstar_gx = bind2.statistic
			bindx_gx = np.linspace(np.min(inl_star_gx[:, 0]), np.max(inl_star_gx[:, 0]), numbins_gx)
			bindy_gx = np.linspace(np.min(inl_star_gx[:, 1]), np.max(inl_star_gx[:, 1]), numbins_gx)
			plt.pcolormesh(bindx_gx, bindy_gx, Nstar_gx, cmap='rainbow', vmin=1e-1, vmax=np.max(Nstar_gx), 
				alpha=1, norm = mpl.colors.LogNorm())
			plt.plot(x0_gx, y0_gx, '+', c = 'b')
			hsc.circles(x0_gx, y0_gx, s = R0_gx, fc = '', ec = 'r')
			plt.title(r'$GadgetX 1.5R_{200} z_{%.3f} Re_{%.3f}$' % ( snap_z_gx, resolution))
			plt.xlim(x0_gx - 1.5*R0_gx, x0_gx + 1.5*R0_gx)
			plt.ylim(y0_gx - 1.5*R0_gx, y0_gx + 1.5*R0_gx)
			plt.xlabel('r[kpc/h]')
			plt.ylabel('r[kpc/h]')
			plt.savefig('GX-z%.3f-distribution.png' % zt, dpi=600)
			plt.show()
			plt.close()
					