import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob
import h5py
import numpy as np
import pandas as pd
import astropy.io.ascii as asc

import find
import changds
import pygadgetreader as pygdr

h = 0.7

C_id = [] # cluster ID
load1 = '/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/'
file1 = glob.glob('/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_*')
C_id = [str(file1[-4:]) for file1 in file1]
CID = np.sort(C_id)

def snap_show():
	load2 = 'NewMDCLUSTER_'
	'''	
	sub_set = ['0110', '0121', '0132', '0182']
	'''
	sub_set = CID
	Nh = len(sub_set)
	for  k in range(Nh):

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

		for kk in range(N):
			zq = kk
			zg = zc[zq]
			Goal_ID = snap[N- zq]
			if zg <= 1:

				try:
					snap_shot_star = pygdr.readsnap(
					'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'pos', 'star')
				except SystemExit:
					print('no star particles now')
					continue
				snap_shot_gas = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'pos', 'gas')
				snap_shot_dm = pygdr.readsnap(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'pos', 'dm')
				snap_N = pygdr.readheader(
				'/mnt/ddnfs/data_users/wgcui/The300/simulation/GadgetX/NewMDCLUSTER_%s/snap_%s' % (sub_set[k], Goal_ID), 'npartTotal')

				halo_list = asc.read(
					'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/GadgetX-NewMDCLUSTER_%s.snap_%s.z%.3f.AHF_halos'
					 % (sub_set[k], sub_set[k], Goal_ID, zg), converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
				Rvir = np.array(halo_list['col12'])
				xhalo = np.array(halo_list['col6'])
				yhalo = np.array(halo_list['col7'])
				zhalo = np.array(halo_list['col8'])

				Halo = np.array(halo_list['col1'])
				sub_halo = np.array(halo_list['col2'])
				goal_halo = HID[trc_id, kk]
				check_id = find.find1d(Halo, goal_halo)

				x0 = xhalo[check_id]
				y0 = yhalo[check_id]
				z0 = zhalo[check_id]
				R0 = Rvir[check_id]

				dr_star = np.sqrt((snap_shot_star[:,0] - x0)**2 + 
					(snap_shot_star[:,1] - y0)**2 + (snap_shot_star[:,2] - z0)**2)
				i_star = dr_star <= R0
				inl_star = snap_shot_star[i_star, :]

				dr_gas = np.sqrt((snap_shot_gas[:,0] - x0)**2 + 
					(snap_shot_gas[:,1] - y0)**2 + (snap_shot_gas[:,2] - z0)**2)
				i_gas = dr_gas <= R0
				inl_gas = snap_shot_gas[i_gas, :]

				dr_dm = np.sqrt((snap_shot_dm[:,0] - x0)**2 + 
					(snap_shot_dm[:,1] - y0)**2 + (snap_shot_dm[:,2] - z0)**2)
				i_dm = dr_dm <= R0
				inl_dm = snap_shot_dm[i_dm, :]

				plt.figure(figsize = (7, 4.5))
				plt.title('Cluster_%s_z%.3f' % (sub_set[k], zg))
				plt.hist2d(inl_gas[:, 0], inl_gas[:,1], bins = [500, 500],
					cmap = 'winter', vmin = 1e-1, vmax = snap_N[0]/10, norm=mpl.colors.LogNorm(), alpha = 0.5)
				plt.hist2d(inl_dm[:, 0], inl_dm[:,1], bins = [500, 500],
					cmap = 'summer', vmin = 1e-1, vmax = snap_N[1]/10, norm=mpl.colors.LogNorm(), alpha = 0.5)
				plt.hist2d(inl_star[:, 0], inl_star[:,1], bins = [500, 500],
					cmap = 'autumn',vmin = 1e-1,vmax = snap_N[-2]/10,norm=mpl.colors.LogNorm(), alpha = 0.5)	
				plt.xlim(x0 - R0, x0 + R0)
				plt.ylim(y0 - R0, y0 + R0)
				#plt.text(x0, y0, s = '%.2f, %.2f'% (x0, y0))
				plt.axis('off')
				plt.xticks([])
				plt.yticks([])
				plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/snap_fig/C_%s_z%.3f.png' % (sub_set[k], zg), dpi = 600)
				plt.close()
			else:
				continue
		print(k)	
	return

def main():
	snap_show()

if __name__ == '__main__':
	main()