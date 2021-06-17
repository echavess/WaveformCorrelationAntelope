import os
import glob
import datetime
import warnings
import numpy as np
import matplotlib as mpl
from obspy.core import read
import matplotlib.pyplot as plt
import cmocean
from mtspec import mt_coherence
import itertools as it
from obspy.io.sac.util import get_sac_reftime
from obspy.signal.filter import lowpass, bandpass, highpass
from obspy.signal.cross_correlation import xcorr_pick_correction
mpl.rcParams['font.sans-serif'] = 'DIN Alternate'
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
warnings.filterwarnings("ignore")


class WCCTD(object):

	def __init__(self, wave_dir, station, trim_time, filter, taus=False):

		self.wave_dir = wave_dir
		self.station = station

		data_path= os.path.join(self.wave_dir, self.station, "*sac")
		sac_files_list = glob.glob(data_path)
		
		Event_id = []
		for sac in sac_files_list:
			event_id = sac.split("/")[2].split("_")[1]
			Event_id.append(event_id)

		st = read(data_path)
		event_id = np.arange(0, len(st))
		print("Total N events in stream: %s" %len(st))
		
		# Combinanción única de eventos para compararœœœœ
		unique = list(it.combinations(np.unique(event_id), 2))
		position_0 = unique[0]
		last_element=(unique[len(unique)-1])
		rows = last_element[0]+1
		cols=last_element[1] + 1

		def get_ptime(trace):
			trace.detrend(type="demean")
			trace.detrend('linear')
			ref_time = get_sac_reftime(trace.stats.sac)
			p_time = ref_time + trace.stats.sac.a
			return trace, p_time,ref_time

		def get_colors(inp, colormap):
			vmin=np.min(inp)
			vmax=np.max(inp)
			norm = plt.Normalize(vmin, vmax)
			return colormap(norm(inp))

		def copying_file(infile, out_directory):
			infile_path = os.path.join(self.wave_dir, self.station, "*%s*" %infile)
			if infile_path is True:
				print("COPYING FILE TO REPS")
				cmd = "cp %s %s" %(infile_path, out_directory)
				os.system(cmd)

		CCs = []
		Taus = []
		Times1 = []
		Times2 = []
		X = []
		Y = []

		file_out = open("correlations_out_%s_%s.dat" %(filter[0], filter[1]), "a")

		for jj, uni in enumerate(unique):
			position = tuple(uni)

			event1_, event2_ = uni[0], uni[1]
			selid1, selid2 = Event_id[uni[0]], Event_id[uni[1]]
			event1, event2 = st[event1_], st[event2_]

			wave1, ptime1, rtime1 = get_ptime(event1)
			wave2, ptime2, rtime2 = get_ptime(event2)

			rtime1 = rtime1.strftime('%b-%dT%H:%M,%Y')
			rtime2 = rtime2.strftime('%b-%dT%H:%M,%Y')

			tau, coeff = xcorr_pick_correction(ptime1, wave1, ptime2, wave2, 
													t_before=trim_time[0],
													t_after=trim_time[1], 
													cc_maxlag=0.8, 
													filter="bandpass",
													filter_options={'freqmin': filterw[0], 'freqmax': filterw[1]},
													plot=False)
			CCs.append(coeff)
			Taus.append(tau)
			Times1.append(rtime1)
			Times2.append(rtime2)
			X.append(event1_)
			Y.append(event2_)

			print("Event pair: %s <---> %s CC: %s, tau: %s, times: %s <---> %s" %(uni[0], uni[1], np.around(coeff, 2), 
											np.around(tau, 2), rtime1, rtime2,))

			if coeff >= 0.95:
				file_out.write("%s %s %s %s %s\n" %(rtime1, rtime2, np.around(coeff, 2), selid1, selid2))
				copying_file(selid1, 'repeaters/')




		file_out.close()

		#### Figure if needed ####
		# CCs = np.asarray(CCs)
		# Taus = np.asarray(Taus)
		# Times1 = list(dict.fromkeys(Times1))
		# Times2 = list(dict.fromkeys(Times2))

		# palette = cmocean.cm.haline
		# Colors = get_colors(CCs, palette)

		# fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8,8))
		# #ax.set_xlim(-0.3, cols-1.5)
		# #ax.set_ylim(0.1, rows+0.5)
		# ax.set_xlabel("Event time")
		# ax.set_ylabel("Event time")
		# ax.yaxis.set_major_locator(MultipleLocator(1))
		# ax.xaxis.set_major_locator(MultipleLocator(1))
		# ax.spines['top'].set_visible(False)
		# ax.spines['right'].set_visible(False)
		# ax.spines['left'].set_visible(False)
		# ax.spines['bottom'].set_visible(False)
		# ax.set_yticks(np.arange(1, rows+1), minor=False)
		# ax.set_xticks(np.arange(0, cols-1), minor=False)
		# ax.set_xticklabels(Times1, fontsize=6, rotation=90)
		# ax.set_yticklabels(Times2, fontsize=6)


		# cc_map = ax.scatter(X, Y, c=CCs, marker="s", s=80, cmap=palette, alpha=0.6, lw=0.7, edgecolor='k', zorder=5)
		# cb = fig.colorbar(cc_map, aspect=30, shrink=0.5)
		# cb.set_label('Correlation coefficient (CC)', rotation=270, labelpad=20, fontsize=11)

		# if taus is True:
		# 	for ii, t in enumerate(Taus):
		# 		ax.text(X[ii]-0.15, Y[ii]+0.15, "%s s" %str(np.around(t, 2)), fontsize=7, color='k')
		# else:
		# 	print("Not including Tau time!")

		#fig.savefig("CCs_map.png", format='png', dpi=700)
		#plt.show()

	
			
wave_dir = 'Waveforms_SAC'
station='JACO'
trim_time=[0.5, 3.5]
filterw=[1, 10]
WCCTD(wave_dir=wave_dir, station=station, 
	trim_time=trim_time, filter=filterw, taus=False)



