"""
trace the BCG galaxy at different snapshot
"""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import h5py
import find
import changds
import numpy as np
import astropy.io.ascii as asc
import pygadgetreader as pygdr
h = 0.7
resolution = 1. # in unit kpc/h
BCG_size = 50 # in unit kpc/h
_id_ = 0
alpha = np.linspace(0, 128, 129)
alpha = np.int0(alpha)
T_scale = len(alpha)
R_in = 50
R_out = 100
def bcg_trace_mu():

	m_bcg = np.zeros(T_scale, dtype = np.float)
	z_mu = np.zeros(T_scale,dtype = np.float)
	with h5py.File('/home/cxkttwl/Scatter/MUSIC/Redshift.h5') as f:
		y0 = f['a']
		com_z = np.array(y0)
	with h5py.File('/home/cxkttwl/Scatter/MUSIC/main_tree.h5') as f:
		y1 = f['a']
		tree_line = np.array(y1)
	L = tree_line.shape[0]
	iv = np.zeros(L,dtype = np.int0)
	for l in range(L):
		u = find.find1d(tree_line[l,:], 0)
		iv[l] = u
	# get the bcg tree line
	z0 = 0.000
	snap_N = pygdr.readheader(
		'/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % '128', 'npartTotal')
	snap_shot_star = pygdr.readsnap(
		'/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % '128', 'pos', 'star')
	main_halo = asc.read(
		'/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z0.000.AHF_halos',
		converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
	goal_z = changds.inv_chidas('%.3f'%z0)
	goalz = find.find1d(com_z, goal_z)
	goal_halo = tree_line[_id_, goalz]
	halo_list = np.array(main_halo['col1'])
	sub_halo = np.array(main_halo['col2'])
	check_id = find.find1d(halo_list, goal_halo)
	Rvir = np.array(main_halo['col12'])
	xhalo = np.array(main_halo['col6'])
	yhalo = np.array(main_halo['col7'])
	zhalo = np.array(main_halo['col8'])
	Mhalo = np.array(main_halo['col4'])
	Mgas = np.array(main_halo['col45'])
	Mstar= np.array(main_halo['col65'])
	# select the total star, gas, dark matter mass
	x0 = xhalo[check_id]
	y0 = yhalo[check_id]
	z0 = zhalo[check_id]
	R0 = Rvir[check_id]

	if snap_N[-2] != 0:
		dstar = np.sqrt((snap_shot_star[:, 0] - x0)**2 + 
		    (snap_shot_star[:, 1] - y0)**2 + (snap_shot_star[:, 2] - z0)**2)
		ids = dstar <= R0
		inlstar = snap_shot_star[ids, :]

	R_range = 150.0
	edge_x = np.array([x0-R_range, x0+R_range])
	edge_y = np.array([y0-R_range, y0+R_range])
	edge_z = np.array([z0-R_range, z0+R_range])
	iv_s = ((inlstar[:,0] <= edge_x[1])&(inlstar[:,0] >= edge_x[0]))&(
		(inlstar[:,1] <= edge_y[1])&(inlstar[:,1] >= edge_y[0]))&(
		(inlstar[:,2] <= edge_z[1])&(inlstar[:,2] >= edge_z[0]))
	inl_star = inlstar[iv_s,:]
	num_bins = np.ceil(R_range*2/resolution)
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

	compare_d= np.sqrt((cen_po_star[:,0]-x0)**2+ (cen_po_star[:,1]-y0)**2+(cen_po_star[:,2]-z0)**2)
	ismin = find.find1d(compare_d,np.min(compare_d))
	cen_x_star = cen_po_star[ismin,0]
	cen_y_star = cen_po_star[ismin,1]
	cen_z_star = cen_po_star[ismin,2]

	ih = sub_halo == halo_list[check_id]
	sub_line = halo_list[ih]
	msat = Mstar[ih]
	sub_x = xhalo[ih]
	sub_y = yhalo[ih]
	sub_z = zhalo[ih]
	sub_r = Rvir[ih]

	dr_sub0 = np.sqrt((sub_x - x0)**2 + (sub_y - y0)**2 + (sub_z - z0)**2)
	isub = dr_sub0 <= R0
	ID_line = sub_line[isub]
	real_sub = msat[isub]
	real_sub_x = sub_x[isub]
	real_sub_y = sub_y[isub]
	real_sub_z = sub_z[isub]
	real_sub_r = sub_r[isub]

	dr_sub1 = np.sqrt((real_sub_x - cen_x_star)**2 + (real_sub_y - cen_y_star)**2 + (real_sub_z - cen_z_star)**2)
	ucen = np.where(dr_sub1 == np.min(dr_sub1))[0]
	trac_ID = ID_line[ucen[0]]
	ix = halo_list == trac_ID
	iy = np.where(ix == True)[0]
	_ID_ = iy[0]

	id_z = com_z[iv[_ID_]]
	print(id_z)
	id_line = tree_line[_ID_]

	for k in range(T_scale):
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
		print('Now redshift is %.3f' % snap_z)
		snap_N = pygdr.readheader(
			'/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
		if snap_z <= id_z:
			z_mu[k] = snap_z
			try:
				snap_shot_star = pygdr.readsnap(
					'/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
				snap_mass_star = pygdr.readsnap(
					'/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'mass','star')
			except SystemExit:
				print('no star particles now')
				break
			main_halo = asc.read(
			        '/mnt/ddnfs/data_users/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos'%snap_z,
			                     converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
			goal_z = changds.inv_chidas('%.3f'%snap_z)
			goalz = find.find1d(com_z,goal_z)
			goal_halo = tree_line[_ID_, goalz]
			halo_list = np.array(main_halo['col1'])
			sub_halo = np.array(main_halo['col2'])
			check_id = find.find1d(halo_list, goal_halo)
			Rvir = np.array(main_halo['col12'])
			xhalo = np.array(main_halo['col6'])
			yhalo = np.array(main_halo['col7'])
			zhalo = np.array(main_halo['col8']) 
			Mhalo = np.array(main_halo['col4'])
			Mgas = np.array(main_halo['col45'])
			Mstar= np.array(main_halo['col65'])
			
	return

def bcg_trace_gx():

	return

def fig():

	return

def main():
	bcg_trace_mu()
	bcg_trace_gx()
	fig()
if __name__ == "__main__":
	main()