import glob
import h5py
import numpy as np
import astropy.io.ascii as asc

C_id = [] # cluster ID
load1 = '/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/'
file1 = glob.glob('/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_*')
C_id = [str(file1[-4:]) for file1 in file1]

def tree_build():

	load2 = 'NewMDCLUSTER_'
	Nh = len(C_id)
	for k in range(Nh):
		sub_file = glob.glob(
			'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/*.AHF_halos' % C_id[k])
		snap_id = [str(sub_file[107:110]) for sub_file in sub_file]
		z = [float(sub_file[112:-10]) for sub_file in sub_file]

		redshift = np.sort(z)
		snap = np.sort(snap_id)
		index_snap = len(snap)-1

		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/z_for_C_%s.h5' % C_id[k],'w') as f:
			f['a'] = np.array(redshift)

		main_value = asc.read(
			'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/GadgetX-NewMDCLUSTER_%s.snap_128.z0.000.AHF_mtree_idx' % (C_id[k], C_id[k]),
			converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
		N = len(main_value['col1'])
		M = len(redshift)
		main_tree = np.zeros((N, M), dtype = np.int64)
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/cluster_%s_tree.h5' % C_id[k],'w') as f:
			f['a'] = np.array(main_tree)

		for q in range(0, len(redshift)-1):
			try:
				subz = redshift[q]
				sub_snap = snap[index_snap - q]
				id_value = asc.read(
				'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/GadgetX-NewMDCLUSTER_%s.snap_%s.z%.3f.AHF_mtree_idx' % (C_id[k], C_id[k], sub_snap, subz),
				converters={'col1':[asc.convert_numpy(np.int64)], 'col2':[asc.convert_numpy(np.int64)]})
				if q == 0 :
					main_tree[:,0] = id_value['col1']
					main_tree[:,1] = id_value['col2']
				else:
					for t in range(len(id_value['col1'])):
						ip_halo = main_tree[:, q] == id_value['col1'][t]
						ipx = ip_halo.tolist()
						ipy = True in ipx
						if ipy == True:
							T_ip = ipx.index(True)
							main_tree[T_ip, q+1] = id_value['col2'][t]
			except FileNotFoundError:
				continue
			print('q = ', q)
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/cluster_%s_tree.h5' % C_id[k]) as f:
			for p in range(len(main_tree)):
				f['a'][p,:] = main_tree[p,:]

	return

def main():
	tree_build()

if __name__ == "__main__":
	main()
