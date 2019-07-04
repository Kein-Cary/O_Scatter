import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.cosmology import Planck15 as Plank
from scipy.interpolate import interp1d as interp

import glob
import h5py
import numpy as np
import astropy.io.ascii as asc

import find
import changds
import pygadgetreader as pygdr

## constant
h = 0.7
C_id = [] # cluster ID
load1 = '/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/'
file1 = glob.glob('/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_*')
C_id = [str(file1[-4:]) for file1 in file1]
CID = np.sort(C_id) ## '0080', '0110' should trace the second one halo at z = 0.

def halo_GX():

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

		zn = np.zeros(N, dtype = np.float)
		mn = np.zeros(N, dtype = np.float)
		ms = np.zeros(N, dtype = np.float)
		mg = np.zeros(N, dtype = np.float)
		print('star select!')
		for q in range(N):
			subz = zc[q]
			subid = snap[N-q]
			if subz <= id_z :
				halo = asc.read(
				'/mnt/ddnfs/data_users/wgcui/The300/catalogues/AHF/GadgetX/NewMDCLUSTER_%s/GadgetX-NewMDCLUSTER_%s.snap_%s.z%.3f.AHF_halos'
				 % (CID[k], CID[k], subid, subz),converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
				print('star read file!')
				Mhalo = np.array(halo['col4'])
				Mgas = np.array(halo['col45'])
				Mstar= np.array(halo['col65'])
				halo_list = np.array(halo['col1'])
				sub_halo = np.array(halo['col2'])

				goal_halo = HID[trc_id, q]
				check_id = find.find1d(halo_list, goal_halo)
				try:
					mn[q] = Mhalo[check_id]
					ms[q] = Mstar[check_id]
					mg[q] = Mgas[check_id]
					zn[q] = subz * 1
				except IndexError:
					continue
			else:
				continue
		ll = ms != 0
		zg = zn[ll]
		Ms = ms[ll]
		Mg = mg[ll]
		Mh = mn[ll]
		record = np.array([zg, Mh, Ms, Mg])
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/mass_reviw_%s_CID.h5' % CID[k], 'w') as f:
			f['a'] = np.array(record)
		with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/mass_reviw_%s_CID.h5' % CID[k]) as f:
			for tt in range(len(record)):
				f['a'][tt,:] = record[tt,:]
		print(k)

	return

def Z_half(z, mh, eta):

	zz = z*1
	mm = mh*1
	A = mm[0] * eta

	fz = interp(mm, zz)
	tz = fz(A)

	dm = np.abs(mm - A)
	idm = np.where(dm == np.min(dm))[0]
	tz0 = zz[idm[0]]
	return tz, tz0

def sample_divide():
	load2 = 'NewMDCLUSTER_'
	Nh = len(CID)
	eta = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
	for tt in range(len(eta)):
		b = eta[tt]
		Z = []
		Mh = []
		Ms = []
		Mg = []
		Z_huf = []
		Z_huf0 = []
		for k in range(Nh):
			with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/mass_reviw_%s_CID.h5' % CID[k]) as f:
				record = np.array(f['a'])
			zc = record[0,:]
			mh = record[1,:]
			ms = record[2,:]
			mg = record[3,:]
			zt, zt0 = Z_half(zc, mh, b)
			Z_huf.append(zt)
			Z_huf0.append(zt0)
			Z.append(zc)
			Mh.append(mh)
			Ms.append(ms)
			Mg.append(mg)

		Zhuf = np.array(Z_huf)
		Zhuf0 = np.array(Z_huf0)
		# sample divided
		Mh0 = [m[0] for m in Mh]
		edg = 0.5 * (np.min(Zhuf) + np.max(Zhuf))
		plt.figure()
		plt.plot(Zhuf, Mh0, 'bo', alpha = 0.5)
		plt.axvline(x = edg, color = 'r', ls = '--')
		plt.xlabel(r'$z_{form} \; of \; halo$')
		plt.ylabel(r'$M^{z=0}_{h} [M_{\odot} / h]$')
		plt.title('Halo mass -- formation time')
		plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/Mass_zform_with%.1f.png' % b, dpi = 300)
		plt.close()

	raise
	return

def fig_out():

	load2 = 'NewMDCLUSTER_'
	Nh = len(CID)
	eta = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
	for tt in range(len(eta)):
		b = eta[tt]
		Z = []
		Mh = []
		Ms = []
		Mg = []
		Z_huf = []
		Z_huf0 = []
		for k in range(Nh):
			with h5py.File('/mnt/ddnfs/data_users/cxkttwl/Scatter/G_x_redshift/tree_h5/mass_reviw_%s_CID.h5' % CID[k]) as f:
				record = np.array(f['a'])
			zc = record[0,:]
			mh = record[1,:]
			ms = record[2,:]
			mg = record[3,:]
			zt, zt0 = Z_half(zc, mh, b)
			Z_huf.append(zt)
			Z_huf0.append(zt0)
			Z.append(zc)
			Mh.append(mh)
			Ms.append(ms)
			Mg.append(mg)

		Zhuf = np.array(Z_huf)
		Zhuf0 = np.array(Z_huf0)
		xzt = np.linspace(0, Nh-1, Nh)
		fig = plt.figure(figsize = (16, 9))
		plt.title('halo mass evolution in Gadget-X')
		ax = plt.subplot(111)
		for p in range(Nh):
			if p == Nh - 1:
				ax.plot(np.log10(1 + Z[p]), Mh[p], ls = '-', color = mpl.cm.rainbow(p/Nh), label = r'$Mh_{cluster %s}$' % CID[p])

				ax1 = ax.twiny()
				xtik = ax.get_xticks()
				Zt = 10**(xtik) - 1
				LBT = Plank.lookback_time(Zt).value
				ax1.set_xticks(xtik)
				ax1.set_xticklabels(["%.2f" % ll for ll in LBT])
				ax1.set_xlim(ax.get_xlim())

			else:
				ax.plot(np.log10(1 + Z[p]), Mh[p], ls = '-', color = mpl.cm.rainbow(p/Nh), label = r'$Mh_{cluster %s}$' % CID[p])
			handles, labels = plt.gca().get_legend_handles_labels()
		ax.legend(loc = 1, fontsize = 10)
		ax.set_yscale('log')
		ax.set_xlabel('$log(1+z)$')
		ax.set_ylabel(r'$Mh[M_{\odot}/h]$')
		ax1.set_yscale('log')
		ax1.set_xlabel(r'$look \; back \; time[Gyr]$')
		ax.tick_params(axis = 'both', which = 'both', direction = 'in')
		ax1.tick_params(axis = 'x', which = 'both', direction = 'in')

		subax = fig.add_axes([0.18, 0.18, 0.42, 0.2])
		subax.set_title(r'$Z_{form} \; with \; %.1f \; M_{h0}$' % b)
		subax.plot(xzt, Zhuf, 'bo', alpha = 0.5)
		subax.plot(xzt, Zhuf0, 'r*', alpha = 0.5)
		subax.set_xlabel('# Cluster')
		subax.set_xticks(xzt)
		subax.set_xticklabels(CID, fontsize = 10)
		subax.set_ylabel(r'$Z_{half}$')

		plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/halo_mass_review_with%.1f.png' % b, dpi = 300)
		plt.close()

	print('halo_assembly_finished!')
	fig = plt.figure(figsize = (16, 9))
	ax = plt.subplot(111)
	plt.title('stellar mass evolution in Gadget-X')
	for p in range(Nh):
		if p == Nh - 1:
			ax.plot(np.log10(1 + Z[p]), Ms[p], ls = '-', color = mpl.cm.rainbow(p/Nh), label = r'$M^{\ast}_{cluster %s}$' % CID[p])

			ax1 = ax.twiny()
			xtik = ax.get_xticks()
			Zt = 10**(xtik) - 1
			LBT = Plank.lookback_time(Zt).value
			ax1.set_xticks(xtik)
			ax1.set_xticklabels(["%.2f" % ll for ll in LBT])
			ax1.set_xlim(ax.get_xlim())

		else:
			ax.plot(np.log10(1 + Z[p]), Ms[p], ls = '-', color = mpl.cm.rainbow(p/Nh), label = r'$M^{\ast}_{cluster %s}$' % CID[p])
		handles, labels = plt.gca().get_legend_handles_labels()
	ax.legend(loc = 1, fontsize = 10)
	ax.set_yscale('log')
	ax.set_xlabel('$log(1+z)$')
	ax.set_ylabel(r'$M_{\ast}[M_{\odot}/h]$')
	ax1.set_yscale('log')
	ax1.set_xlabel(r'$look \; back \; time[Gyr]$')
	ax.tick_params(axis = 'both', which = 'both', direction = 'in')
	ax1.tick_params(axis = 'x', which = 'both', direction = 'in')	

	plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/stellar_mass_review.png', dpi = 300)
	plt.close()

	fig = plt.figure(figsize = (16, 9))
	ax = plt.subplot(111)
	plt.title('gass mass evolution in Gadget-X')
	for p in range(Nh):
		if p == Nh - 1:
			ax.plot(np.log10(1 + Z[p]), Mg[p], ls = '-', color = mpl.cm.rainbow(p/Nh), label = r'$M^{gas}_{cluster %s}$' % CID[p])

			ax1 = ax.twiny()
			xtik = ax.get_xticks()
			Zt = 10**(xtik) - 1
			LBT = Plank.lookback_time(Zt).value
			ax1.set_xticks(xtik)
			ax1.set_xticklabels(["%.2f" % ll for ll in LBT])
			ax1.set_xlim(ax.get_xlim())

		else:
			ax.plot(np.log10(1 + Z[p]), Mg[p], ls = '-', color = mpl.cm.rainbow(p/Nh), label = r'$M^{gas}_{cluster %s}$' % CID[p])
		handles, labels = plt.gca().get_legend_handles_labels()
	ax.legend(loc = 1, fontsize = 10)
	ax.set_yscale('log')
	ax.set_xlabel('$log(1+z)$')
	ax.set_ylabel(r'$M_{gas}[M_{\odot}/h]$')
	ax1.set_yscale('log')
	ax1.set_xlabel(r'$look \; back \; time[Gyr]$')
	ax.tick_params(axis = 'both', which = 'both', direction = 'in')
	ax1.tick_params(axis = 'x', which = 'both', direction = 'in')

	plt.savefig('/mnt/ddnfs/data_users/cxkttwl/Scatter/snap/gas_mass_review.png', dpi = 300)
	plt.close()

	return

def main():
	#halo_GX()
	sample_divide()
	#fig_out()

if __name__ == "__main__":
	main()