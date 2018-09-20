# this file try to read the snapshot data to find the particle information
import matplotlib as mpl
mpl.use('Agg')
#from pygadgetreader import *
import matplotlib.pyplot as plt
#from handy import scatter as hsc
import numpy as np
import h5py
import pygadgetreader as pygdr
#import pandas as pd
import astropy.io.ascii as asc
import find
import changds
def fig_snap(ip,resolution,viewp):
    """
    ip : the goal halo number of halo at z==0
    resolution : the resolution of graph
    viewp : the view direction,must be integer
    # view control the projection of distribution：1--x-y panel,2--x-z panel,3--y-z panel
    """
    _id_ = ip
    resolution = resolution
    alpha = np.linspace(0, 128, 129)
    alpha = np.int0(alpha)
    # section2:read out the position data of MUSIC
    plane_id = np.int0(viewp)
    view = np.array([1,2,3])
    view_id = view[plane_id-1]
    #read the meger tree to get the star time of each halo
    with h5py.File('/home/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/Redshift.h5') as f:
        y0 = f['a']
        com_z = np.array(y0)
    with h5py.File('/home/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/main_tree.h5') as f:
        y1 = f['a']
        tree_line = np.array(y1)
    with h5py.File('/home/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/Redshift_GX.h5') as f:
        y2 = f['a']
        com_z_gx = np.array(y2)
    with h5py.File('/home/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/main_tree_GX.h5') as f:
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
        #print('Now redshift is %.3f'%snap_z)
        snap_N = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
        # print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
        #      %(snap_N[0],snap_N[1],snap_N[2],snap_N[3],snap_N[4],snap_N[5]))
        snap_name = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'header')
        # for type
        if snap_z <= id_z:
            snap_shot_gas = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_shot_DM = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
            try:
                snap_shot_star = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
            except SystemExit:
                print('no star particles now')
            # for respective position
            snap_shot_bulge = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bulge')
            snap_shot_disk = pygdr.readsnap(
                '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'disk')
            try:
                snap_shot_bndry = pygdr.readsnap(
                    '/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bndry')
            except SystemExit:
                print('no boundary particles now')
            # for density profile of gas
            #snap_dens = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetMUSIC/NewMDCLUSTER_0001/snap_%s'%No_snap,'rho','gas')
            '''
            # the total mass distribution of snap_shot_128
            plt.figure()
            if snap_N[-2] != 0:
                plt.hist2d(snap_shot_star[:, 0], snap_shot_star[:, 1], bins=[1000, 1000],
                           cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
                plt.axis('off')
                plt.xticks([])
                plt.yticks([])
                plt.title('MUSIC tot star z_%.3f' % snap_z)
                plt.savefig('snap_fig/MUSIC star distribution z_%.0f.png' % dd, dpi=600)
                plt.show()
                plt.close()
            plt.figure()
            plt.hist2d(snap_shot_gas[:, 0], snap_shot_gas[:, 1], bins=[1000, 1000],
                       cmap='rainbow', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
            plt.colorbar(label=r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('MUSIC tot gas z_%.3f' % snap_z)
            plt.savefig('snap_fig/MUSIC gas distribution z_%.0f.png' % dd, dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            plt.hist2d(snap_shot_DM[:, 0], snap_shot_DM[:, 1], bins=[1000, 1000],
                       cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
            plt.colorbar(label=r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('MUSIC tot DM z_%.3f' % snap_z)
            plt.savefig('snap_fig/MUSIC DM distribution z_%.0f.png' % dd, dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            if snap_N[-1] !=0:
                plt.hist2d(snap_shot_bndry[:,0],snap_shot_bndry[:,1],bins = [1000,1000],
                           cmap = 'GnBu',vmin = 1e-1,vmax = snap_N[-1]/100,norm=mpl.colors.LogNorm())
                plt.colorbar(label = r'$log[Number density]$')
                plt.axis('off')
                plt.xticks([])
                plt.yticks([])
                plt.title('MUSIC tot bnd z_%.3f' % snap_z)
                plt.savefig('snap_fig/MUSIC boundary distribution z_%.0f.png'% dd,dpi=600)
                plt.show()
                plt.close()
            plt.figure()
            plt.hist2d(snap_shot_disk[:,0],snap_shot_disk[:,1],bins = [1000,1000],
                       cmap = 'cool',vmin = 1e-1,vmax = snap_N[2]/100,norm=mpl.colors.LogNorm())
            plt.colorbar(label = r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('MUSIC tot disk z_%.3f'% snap_z)
            plt.savefig('snap_fig/MUSIC disk distribution z_%.0f.png'% dd,dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            plt.hist2d(snap_shot_bulge[:,0],snap_shot_bulge[:,1],bins = [1000,1000],
                       cmap = 'summer',vmin = 1e-1,vmax = snap_N[3]/100,norm=mpl.colors.LogNorm())
            plt.colorbar(label = r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('MUSIC tot bul z_%.3f'% snap_z)
            plt.savefig('snap_fig/MUSIC bulge distribution z_%.0f.png'% dd,dpi=600)
            plt.show()
            plt.close()
            '''
            # next,try to find the particles belong to halo 128000000000001(in MUSIC simulation)
            main_halo = asc.read('/home/cxkttwl/Scatter_data_read/MUSIC_reshift/NewMDCLUSTER_0001/GadgetMUSIC-NewMDCLUSTER_0001.z%.3f.AHF_halos' % snap_z,
                                 converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z = changds.inv_chidas('%.3f'%snap_z)
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
            ig = dgas <= R0
            inl_gas = snap_shot_gas[ig, :]
            dDM = np.sqrt((snap_shot_DM[:, 0]-x0)**2+(snap_shot_DM[:, 1]-y0)**2 +
                          (snap_shot_DM[:, 2]-z0)**2)
            iD = dDM <= R0
            inl_DM = snap_shot_DM[iD, :]
            ddisk = np.sqrt((snap_shot_disk[:, 0]-x0)**2+(snap_shot_disk[:, 1]-y0)**2 +
                            (snap_shot_disk[:, 2]-z0)**2)
            idisk = ddisk <= R0
            inl_disk = snap_shot_disk[idisk, :]
            dbulge = np.sqrt((snap_shot_bulge[:, 0]-x0)**2+(snap_shot_bulge[:, 1]-y0)**2 +
                             (snap_shot_bulge[:, 2]-z0)**2)
            ibu = dbulge <= R0
            inl_bulge = snap_shot_bulge[ibu, :]
            if snap_N[-2] != 0:
                dstar = np.sqrt((snap_shot_star[:, 0]-x0)**2+(snap_shot_star[:, 1]-y0)**2 +
                                (snap_shot_star[:, 2]-z0)**2)
                ids = dstar <= R0
                inl_star = snap_shot_star[ids, :]
            if snap_N[-1] != 0:
                dbnd = np.sqrt((snap_shot_bndry[:, 0]-x0)**2+(snap_shot_bndry[:, 1]-y0)**2 +
                               (snap_shot_bndry[:, 2]-z0)**2)
                ibnd = dbnd <= R0
                inl_bndry = snap_shot_bndry[ibnd, :]
            numbins = np.ceil(2*R0/resolution)
            plt.figure()
            if view_id == 1:
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 1], bins=[numbins, numbins],
                       cmap='rainbow', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 2:
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 2], bins=[numbins, numbins],
                       cmap='rainbow', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 3:
                plt.hist2d(inl_gas[:, 1], inl_gas[:, 2], bins=[numbins, numbins],
                       cmap='rainbow', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            plt.xlim(x0-R0, x0+R0)
            plt.ylim(y0-R0, y0+R0)
            # plt.axis('off')
            # plt.xticks([])
            # plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.xticks([x0-R0, x0, x0+R0], size=10)
            plt.yticks([y0, y0+R0], rotation=90, size=10)
            plt.title(r'$MUSIC_{%.0f}-inl-gas-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z, resolution))
            plt.savefig('snap_fig/MUSIC gas distribution No_%.0f_z_%.0f.png' %(_id_, dd), dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            if view_id == 1:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 1], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 2:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 3:
                plt.hist2d(inl_DM[:, 1], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            plt.xlim(x0-R0, x0+R0)
            plt.ylim(y0-R0, y0+R0)
            # plt.axis('off')
            # plt.xticks([])
            # plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.xticks([x0-R0, x0, x0+R0], size =10)
            plt.yticks([y0, y0+R0], rotation=90, size =10)
            plt.title(r'$MUSIC_{%.0f}-inl-DM-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z, resolution))
            plt.savefig('snap_fig/MUSIC DM distribution No_%.0f_z_%.0f.png' %(_id_, dd), dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            if snap_N[-2] != 0:
                if view_id == 1:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 1], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                if view_id == 2:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                if view_id == 3:
                    plt.hist2d(inl_star[:, 1], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                plt.scatter(x0, y0, s=50, c='k', marker='+')
                plt.xlim(x0-R0, x0+R0)
                plt.ylim(y0-R0, y0+R0)
                # plt.axis('off')
                # plt.xticks([])
                # plt.yticks([])
                plt.xlabel(r'$r-kpc/h$')
                plt.ylabel(r'$r-kpc/h$')
                plt.xticks([x0-R0, x0, x0+R0], size=10)
                plt.yticks([y0, y0+R0], rotation=90, size=10)
                plt.title(r'$MUSIC_{%.0f}-inl-star-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z, resolution))
                plt.savefig('snap_fig/MUSIC star distribution No_%.0f_z_%.0f.png' %(_id_, dd), dpi=600)
                plt.show()
                plt.close()
            # multi-distribution
            plt.figure(figsize=(8, 8))
            if view_id == 1:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 1], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 1], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 1], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM[:, 1], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 1], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 1], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            plt.xlim(x0-R0, x0+R0)
            plt.ylim(y0-R0, y0+R0)
            # plt.axis('off')
            # plt.xticks([])
            # plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetMUSIC_{%.0f}-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z, resolution))
            plt.savefig('snap_fig/MUSIC Multi-distribution No_%.0f_z_%.0f.png' %(_id_, dd), dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            if view_id == 1:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 1], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 1], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 1], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM[:, 1], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 1], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 1], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            plt.scatter(x0, y0, s=50, c='k', marker='x',)
            plt.xlim(x0-R0*3.5/22.0, x0+R0*3.5/22.0)
            plt.ylim(y0-R0*3.0/14.2, y0+R0*3.0/14.2)
            plt.xticks([x0-R0*3.5/22.0, x0, x0+R0*3.5/22.0], size=10)
            plt.yticks([y0, y0+R0*3.0/14.2], rotation=90, size=10)
            # plt.axis('off')
            # plt.xticks([])
            # plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetMUSIC_{%.0f}-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z, resolution))
            plt.savefig('snap_fig/MUSIC Multi-center No_%.0f_z_%.0f.png' % (_id_, dd), dpi=600)
            plt.show()
            plt.close()
        # section3:read out the position data of GadgetX
        # for GX simulation
        id_z_gx = com_z_gx[iv_gx[_id_]]
        snap_z_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'redshift')
        snap_z_gx = np.abs(snap_z_gx)
        #print('Now redshift is %.3f'%snap_z_gx)
        snap_N_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'npartTotal')
        # print('Now particles:\n gas %.0f\n dark matter %.0f\n disk %.0f\n bulge %.0f\n star %.0f\n boundary %.0f'
        #      %(snap_N_gx[0],snap_N_gx[1],snap_N_gx[2],snap_N_gx[3],snap_N_gx[4],snap_N_gx[5]))
        snap_name_gx = pygdr.readheader('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'header')
        # for type
        if snap_z_gx <= id_z_gx:
            snap_shot_gas_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'gas')
            snap_shot_DM_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'dm')
            try:
                snap_shot_star_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'star')
            except SystemExit:
                print('no star particles now')
            # for respective position
            snap_shot_bulge_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bulge')
            snap_shot_disk_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'disk')
            try:
                snap_shot_bndry_gx = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s' % No_snap, 'pos', 'bndry')
            except SystemExit:
                print('no boundary particles now')
            # for density profile of gas
            #snap_dens = pygdr.readsnap('/mnt/ddnfs/data_users/wgcui/The300/GadgetX/NewMDCLUSTER_0001/snap_%s'%No_snap,'rho','gas')
            '''
            # the total mass distribution of snap_shot_128
            plt.figure()
            if snap_N_gx[-2] != 0:
                plt.hist2d(snap_shot_star_gx[:, 0], snap_shot_star_gx[:, 1], bins=[1000, 1000],
                           cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
                plt.axis('off')
                plt.xticks([])
                plt.yticks([])
                plt.title('GX tot star z_%.3f' % snap_z_gx)
                plt.savefig('snap_fig/GX star distribution z_%.0f.png' % dd, dpi = 600)
                plt.show()
                plt.close()
            plt.figure()
            plt.hist2d(snap_shot_gas_gx[:, 0], snap_shot_gas_gx[:, 1], bins=[1000, 1000],
                       cmap='rainbow', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
            plt.colorbar(label=r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('GX tot gas z_%.3f' % snap_z_gx)
            plt.savefig('snap_fig/GX gas distribution z_%.0f.png' % dd, dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            plt.hist2d(snap_shot_DM_gx[:, 0], snap_shot_DM_gx[:, 1], bins=[1000, 1000],
                       cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
            plt.colorbar(label=r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('GX tot DM z_%.3f' % snap_z_gx)
            plt.savefig('snap_fig/GX DM distribution z_%.0f.png' % dd, dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            if snap_N_gx[-1] !=0:
                plt.hist2d(snap_shot_bndry_gx[:,0],snap_shot_bndry_gx[:,1],bins = [500,500],
                           cmap = 'GnBu',vmin = 1e-1,vmax = snap_N_gx[-1]/100,norm=mpl.colors.LogNorm())
                plt.colorbar(label = r'$log[Number density]$')
                plt.axis('off')
                plt.xticks([])
                plt.yticks([])
                plt.title('GX tot bnd z_%.3f' %snap_z_gx)
                plt.savefig('snap_fig/GX boundary distribution z_%.0f.png'%dd,dpi=600)
                plt.show()
                plt.close()
            plt.figure()
            plt.hist2d(snap_shot_disk_gx[:,0],snap_shot_disk_gx[:,1],bins = [500,500],
                       cmap = 'cool',vmin = 1e-1,vmax = snap_N_gx[2]/100,norm=mpl.colors.LogNorm())
            plt.colorbar(label = r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('GX tot disk z_%.3f'% snap_z_gx)
            plt.savefig('snap_fig/GX disk distribution z_%.0f.png'%dd,dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            plt.hist2d(snap_shot_bulge_gx[:,0],snap_shot_bulge_gx[:,1],bins = [500,500],
                       cmap = 'summer',vmin = 1e-1,vmax = snap_N_gx[3]/100,norm=mpl.colors.LogNorm())
            plt.colorbar(label = r'$log[Number density]$')
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.title('GX tot bulg z_%.3f'% snap_z_gx)
            plt.savefig('snap_fig/GX bulge distribution z_%.0f.png'%dd,dpi=600)
            plt.show()
            plt.close()
            '''
            # next,try to find the particles belong to halo 128000000000001(in GX simulation)
            main_halo_gx = asc.read('/home/cxkttwl/Scatter_data_read/G_x_redshift/NewMDCLUSTER_0001/GadgetX-NewMDCLUSTER_0001.z%.3f.AHF_halos' % snap_z_gx,
                                    converters={'col1': [asc.convert_numpy(np.int64)], 'col2': [asc.convert_numpy(np.int64)]})
            goal_z_gx = changds.inv_chidas('%.3f'%snap_z_gx)
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
            ig_gx = dgas_gx <= R0_gx
            inl_gas_gx = snap_shot_gas_gx[ig_gx, :]
    
            dDM_gx = np.sqrt((snap_shot_DM_gx[:, 0]-x0_gx)**2+(snap_shot_DM_gx[:, 1]-y0_gx)**2 +
                             (snap_shot_DM_gx[:, 2]-z0_gx)**2)
            iD_gx = dDM_gx <= R0_gx
            inl_DM_gx = snap_shot_DM_gx[iD_gx, :]
            ddisk_gx = np.sqrt((snap_shot_disk_gx[:, 0]-x0_gx)**2+(snap_shot_disk_gx[:, 1]-y0_gx)**2 +
                               (snap_shot_disk_gx[:, 2]-z0_gx)**2)
            idisk_gx = ddisk_gx <= R0_gx
            inl_disk_gx = snap_shot_disk_gx[idisk_gx, :]
            dbulge_gx = np.sqrt((snap_shot_bulge_gx[:, 0]-x0_gx)**2+(snap_shot_bulge_gx[:, 1]-y0_gx)**2 +
                                (snap_shot_bulge_gx[:, 2]-z0_gx)**2)
            ibu_gx = dbulge_gx <= R0_gx
            inl_bulge_gx = snap_shot_bulge_gx[ibu_gx, :]
            if snap_N_gx[-2] != 0:
                dstar_gx = np.sqrt((snap_shot_star_gx[:, 0]-x0_gx)**2+(snap_shot_star_gx[:, 1]-y0_gx)**2 +
                                   (snap_shot_star_gx[:, 2]-z0_gx)**2)
                ids_gx = dstar_gx <= R0_gx
                inl_star_gx = snap_shot_star_gx[ids_gx, :]
            if snap_N_gx[-1] != 0:
                dbnd_gx = np.sqrt((snap_shot_bndry_gx[:, 0]-x0_gx)**2+(snap_shot_bndry_gx[:, 1]-y0_gx)**2 +
                                  (snap_shot_bndry_gx[:, 2]-z0_gx)**2)
                ibnd_gx = dbnd_gx <= R0_gx
                inl_bndry_gx = snap_shot_bndry_gx[ibnd_gx, :]
            numbins_gx = np.ceil(2*R0_gx/resolution)
            plt.figure()
            if view_id == 1:
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='rainbow', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 2:
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='rainbow', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 3:
                plt.hist2d(inl_gas_gx[:, 1], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='rainbow', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            plt.xlim(x0_gx-R0_gx, x0_gx+R0_gx)
            plt.ylim(y0_gx-R0_gx, y0_gx+R0_gx)
            plt.xticks([x0_gx-R0_gx, x0_gx, x0_gx+R0_gx], size=10)
            plt.yticks([y0_gx, y0_gx+R0_gx], rotation=90, size=10)
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GX_{%.0f}-inl-gas-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z_gx, resolution))
            plt.savefig('snap_fig/GX gas distribution No_%.0f_z_%.0f.png' % (_id_, dd), dpi = 600)
            plt.show()
            plt.close()
            plt.figure()
            if view_id == 1:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 2:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            if view_id == 3:
                plt.hist2d(inl_DM_gx[:, 1], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.colorbar(label=r'$log[Number density]$')
            plt.xlim(x0_gx-R0_gx, x0_gx+R0_gx)
            plt.ylim(y0_gx-R0_gx, y0_gx+R0_gx)
            plt.xticks([x0_gx-R0_gx, x0_gx, x0_gx+R0_gx], size=10)
            plt.yticks([y0_gx, y0_gx+R0_gx], rotation=90, size=10)
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GX_{%.0f}-inl-DM-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z_gx, resolution))
            plt.savefig('snap_fig/GX DM distribution No_%.0f_z_%.0f.png' % (_id_, dd), dpi = 600)
            plt.show()
            plt.close()
            plt.figure()
            if snap_N[-2] != 0:
                if view_id == 1:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 1], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                if view_id == 2:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                if view_id == 3:
                    plt.hist2d(inl_star_gx[:, 1], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
                    plt.colorbar(label=r'$log[Number density]$')
                plt.scatter(x0_gx, y0_gx, s=50, c='k', marker='+')
                plt.xlim(x0_gx-R0_gx, x0_gx+R0_gx)
                plt.ylim(y0_gx-R0_gx, y0_gx+R0_gx)
                plt.xticks([x0_gx-R0_gx, x0_gx, x0_gx+R0_gx], size=10)
                plt.yticks([y0_gx, y0_gx+R0_gx], rotation=90, size=10)
                plt.xlabel(r'$r-kpc/h$')
                plt.ylabel(r'$r-kpc/h$')
                plt.title(r'$GX_{%.0f}-inl-star-z_{%.3f} Re_{%.3f}$' % (plane_id,snap_z_gx, resolution))
                plt.savefig('snap_fig/GX star distribution No_%.0f_z_%.0f.png' % (_id_, dd), dpi = 600)
                plt.show()
                plt.close()
            # multi-distribution
            plt.figure(figsize=(8, 8))
            if view_id == 1:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 1], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM_gx[:, 1], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 1], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 1], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            plt.xlim(x0_gx-R0_gx, x0_gx+R0_gx)
            plt.ylim(y0_gx-R0_gx, y0_gx+R0_gx)
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetX_{%.0f}-z_{%.3f} Re_{%.3f}$' % (plane_id, snap_z_gx, resolution))
            plt.savefig('snap_fig/GX Multi-distribution No_%.0f_z_%.0f.png' % (_id_, dd), dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            if view_id == 1:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 1], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM_gx[:, 1], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 1], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 1], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            plt.scatter(x0_gx, y0_gx, s=50, c='k', marker='x',)
            plt.xlim(x0_gx-R0_gx*3.5/22.0, x0_gx+R0_gx*3.5/22.0)
            plt.ylim(y0_gx-R0_gx*3.0/14.2, y0_gx+R0_gx*3.0/14.2)
            plt.xticks([x0_gx-R0_gx*3.5/22.0, x0_gx,
                        x0_gx+R0_gx*3.5/22.0], size=10)
            plt.yticks([y0_gx, y0_gx+R0_gx*3.0/14.2], rotation=90, size=10)
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetX_{%.0f}-z_{%.3f} Re_{%.3f}$' % (plane_id, snap_z_gx, resolution))
            plt.savefig('snap_fig/GX Multi-center No_%.0f_z_%.0f.png' % (_id_, dd), dpi=600)
            plt.show()
            plt.close()
        if snap_z <= id_z and snap_z_gx <= id_z_gx:
            # comparation from x-z panel
            plt.figure(figsize=(8, 8))
            plt.subplot(1, 2, 1)
            if view_id == 1:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 1], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 1], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 1], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM[:, 1], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 1], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 1], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            plt.xlim(x0-R0, x0+R0)
            plt.ylim(y0-R0, y0+R0)
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetMUSIC_{%.0f}-z_{%.3f}-Re_{%.3f}$' % (plane_id,snap_z, resolution))
            plt.subplot(1, 2, 2)
            if view_id == 1:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 1], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM[:, 1], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 1], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 1], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            plt.xlim(x0_gx-R0_gx, x0_gx+R0_gx)
            plt.ylim(y0_gx-R0_gx, y0_gx+R0_gx)
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetX_{%.0f}-z_{%.3f}-Re_{%.3f}$' % (plane_id,snap_z_gx, resolution))
            plt.tight_layout()
            plt.savefig('snap_fig/Image halo_%.0f particle snap_%.0f.png' %(_id_,dd), dpi=600)
            plt.show()
            plt.close()
            plt.figure()
            plt.subplot(1, 2, 1)
            if view_id == 1:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 1], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 1], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 1], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM[:, 0], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 0], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 0], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM[:, 1], inl_DM[:, 2], bins=[numbins, numbins],
                           cmap='viridis', vmin=1e-1, vmax=snap_N[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas[:, 1], inl_gas[:, 2], bins=[numbins, numbins],
                           cmap='cool', vmin=1e-1, vmax=snap_N[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star[:, 1], inl_star[:, 2], bins=[numbins, numbins],
                               cmap='plasma', vmin=1e-1, vmax=snap_N[-2]/100, norm=mpl.colors.LogNorm())
            plt.scatter(x0, y0, s=50, c='k', marker='x',)
            plt.xlim(x0-R0*3.5/22.0, x0+R0*3.5/22.0)
            plt.ylim(y0-R0*3.0/14.2, y0+R0*3.0/14.2)
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetMUSIC_{%.0f}-z_{%.3f}-Re_{%.3f}$' % (plane_id,snap_z, resolution))
            plt.subplot(1, 2, 2)
            if view_id == 1:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 1], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 1], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 2:
                plt.hist2d(inl_DM_gx[:, 0], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 0], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 0], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            if view_id == 3:
                plt.hist2d(inl_DM_gx[:, 1], inl_DM_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='viridis', vmin=1e-1, vmax=snap_N_gx[1]/100, norm=mpl.colors.LogNorm())
                plt.hist2d(inl_gas_gx[:, 1], inl_gas_gx[:, 2], bins=[numbins_gx, numbins_gx],
                           cmap='cool', vmin=1e-1, vmax=snap_N_gx[0]/100, norm=mpl.colors.LogNorm())
                if snap_N[-2] != 0:
                    plt.hist2d(inl_star_gx[:, 1], inl_star_gx[:, 2], bins=[numbins_gx, numbins_gx],
                               cmap='plasma', vmin=1e-1, vmax=snap_N_gx[-2]/100, norm=mpl.colors.LogNorm())
            plt.scatter(x0_gx, y0_gx, s=50, c='k', marker='x',)
            plt.xlim(x0_gx-R0_gx*3.5/22.0, x0_gx+R0_gx*3.5/22.0)
            plt.ylim(y0_gx-R0_gx*3.0/14.2, y0_gx+R0_gx*3.0/14.2)
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.xlabel(r'$r-kpc/h$')
            plt.ylabel(r'$r-kpc/h$')
            plt.title(r'$GadgetX_{%.0f}-z_{%.3f}-Re_{%.3f}$' % (plane_id,snap_z_gx, resolution))
            plt.tight_layout()
            plt.savefig('snap_fig/Image halo_%.0f center snap_%.0f.png' % (_id_,dd), dpi=600)
            plt.show()
            plt.close()
