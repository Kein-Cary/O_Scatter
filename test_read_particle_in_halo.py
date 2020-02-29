import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob
import h5py
import struct
import numpy as np
import astropy.io.ascii as asc
import pygadgetreader as pygdr

def read_particle():
	load1 = 'D:/mask/snapshot/GX/snap_128'
	comp = pygdr.readheader(load1, 'npartTotal')
	f = open(load1, 'rb')
	for Nx in range(10):
		for Ny in range(10):
			data = f.read(4)
			head = struct.unpack("i", data)[0]
			print(head)
	raise
	return

def main():
	read_particle()

if __name__ == "__main__":
	main()