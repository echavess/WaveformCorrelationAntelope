''' 
Waveform_relocation code for Double difference earthquake relocation
Author: Dr. Esteban J. Chaves
Volcanological and Seismological Observatory of Costa Rica
OVSICORI-UNA
2019
'''
import glob
import os, shutil
import numpy as np
import datetime
import itertools as it
import pandas as pd
from obspy.core import read
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from obspy import UTCDateTime
from obspy.io.sac.util import get_sac_reftime
from obspy.signal.filter import lowpass, bandpass, highpass
from obspy.signal.cross_correlation import xcorr_pick_correction
from obspy.core.util.attribdict import AttribDict
from eqcorrscan.utils.clustering import cluster, distance_matrix
import warnings
warnings.simplefilter("ignore")

# Colormaps
from matplotlib import cm as cm

# Clustering/Hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

import seaborn as sns
sns.set()

from .prep_event_list import _prep_in_evtlist
from .unify_xcorr_tables import unify_xcorr_tables

class Waveform_relocation(object):

	def __init__(self, Miniseed_dir, Seismic_Catalogs, min_mag, max_mag):

		print("----- Waveform relocation codes developed @ OVSICORI-UNA by Dr. Esteban J. Chaves -----")

		self.Miniseed_dir = Miniseed_dir
		self.event_locations = Seismic_Catalogs[0]
		self.event_phases = Seismic_Catalogs[1]

		self.min_mag = min_mag
		self.max_mag = max_mag

		if not os.path.exists(self.Miniseed_dir):
			_msg0 = self.Miniseed_dir + " does not exists, please check path"
			raise IOError(_msg0)

		if not os.path.isfile(self.event_locations):
			_msg00 = "The seismic Catalog declared is not in current directory, please check name or path"
			raise IOError(_msg00)


	def get_templates_from_catalog(self, station, stime, endtime, save_waves):

		Day, Time, Latitude, Longitude, Depth, Magnitude, evid = np.genfromtxt(self.event_locations, unpack=True, 
																			   usecols=(0, 1, 2, 3, 4, 5, 9,), 
																			   dtype='str')
											
		Magnitude = Magnitude.astype(np.float)
		Latitude = Latitude.astype(np.float)
		Longitude = Longitude.astype(np.float)
		Depth = Depth.astype(np.float)
		evid = evid.astype(np.int)

		stime = stime
		etime = endtime

		Time_vec = []
		for i, day in enumerate(Day):
			temp = day + "T" + Time[i]
			Time_vec.append(datetime.datetime.strptime(temp, "%Y%m%dT%H%M%S%f"))
		
		Time_vec = np.asarray(Time_vec)
		

		sel_mag = np.where(np.logical_and(Magnitude >= self.min_mag, Magnitude <= self.max_mag))
		Magnitude = Magnitude[sel_mag]
		Latitude = Latitude[sel_mag]
		Longitude = Longitude[sel_mag]
		Depth = Depth[sel_mag]
		Time_vec = Time_vec[sel_mag]
		evid = evid[sel_mag]

		print("Amount of events detected in this magnitude range: %s" %len(evid), evid)

		p_times = []
		
		for iids in evid:
			index = "#  %s" %iids
			print("extracting phase time for event %s, station: %s" %(index, station))

			stations_used = []

			with open(self.event_phases) as f:
				for line in f:
					if line.startswith(index):
						for line in f:
							data = line.split(' ')
							data = list(filter(None, data))
							if len(data) == 4:
								if data[3] == 'P\n':
									stations_used.append(data[0])
							
							else:
								break

			if not station in stations_used:
				print("Station was not used for this event...")
				p_times.append('nan')
			
			else:
				with open(self.event_phases) as f:
					for line in f:
						if line.startswith(index):
							for line in f:
								data = line.split(' ')
								data = list(filter(None, data))
								if len(data) == 4:
									if data[0] ==  station:
										if data[3] == 'P\n':
											p_times.append(float(data[1])+stime)
								
								else: 
									break
				
		p_times = np.asarray(p_times)
													
		out_dir = "Waveforms_%s" %self.min_mag + "_%s" %self.max_mag
		stat_dir = os.path.join(out_dir, station)

		if not os.path.exists(out_dir):
			os.mkdir(out_dir)
		
		else:
			print(out_dir + " it already exists!!")

		if not os.path.exists(stat_dir):
			os.mkdir(stat_dir)

		else:
			print(stat_dir + " it already exists!!")
		
		text_file_out = "Waveforms_%s" %self.min_mag + "_%s" %self.max_mag + '.dat'
		text_out_path ="%s/%s" %(stat_dir, text_file_out)
		outfile = open(text_out_path,'a')


		for ii, t in enumerate(Time_vec): 
			jul = t.strftime('%j')
			year = t.year
			jul_ = float(t.strftime('%j'))
			
			dt = UTCDateTime(t)
			#station_to_read = self.Miniseed_dir + '/OV.%s.HHZ..%s.%s.msd' %(station, year, jul)
			station_to_read = self.Miniseed_dir + '/i4.%s' %station + '.HHZ.%s'%year + jul +'_0+'

			if os.path.exists(station_to_read):

				st = read(station_to_read, starttime=dt-stime, endtime=dt+etime)
				for i, tr in enumerate(st):
					staname = tr.stats.station
					chan = tr.stats.channel
					time= tr.stats.starttime
					year = str(time.year)
					jday = str(time.julday)

					tr.stats.sac = AttribDict()
					tr.stats.sac.o = 0
					tr.stats.sac.evdp = 0 
					tr.stats.sac.lovrok = 1
					tr.stats.sac.lcalda = 1
					tr.stats.sac.evla = Latitude[ii]
					tr.stats.sac.evlo = Longitude[ii]
					tr.stats.sac.evdp = Depth[ii]
					tr.stats.sac.mag = Magnitude[ii]
					tr.stats.sac.nzyear = year
					tr.stats.sac.nzjday = jday
					# tr.stats.sac.nzhour = Hour
					# tr.stats.sac.nzmin = Min
					# tr.stats.sac.nzsec = Sec
					# tr.stats.sac.nzmsec = 0
					tr.stats.sac.knetwk = 'OV'

					ppick = p_times[ii]
					if ppick != 'nan':
						tr.stats.sac.a = p_times[ii]

					# Filename and new output location
					file_name = staname+"."+chan+"."+ year +"."+str(time)+ "_%s_" %evid[ii] + ".sac"
					new_loc = os.path.join(out_dir, station, file_name)
					outfile.write("%s\n" %(evid[ii]))
					
					# Wirtting the waveforms
					print("Working on file = %s" %file_name)
					tr.write(new_loc, format="SAC")


			else:
				print("%s doesn't exists in directory..." %station_to_read)
				continue
		
		outfile.close()

		path = os.path.join(stat_dir, '*sac')
		stream_files = glob.glob(path)
		stream_list = [(read(stream_file), i) for i, stream_file in enumerate(stream_files)]
		groups = cluster(template_list=stream_list, show=False, corr_thresh=0.3, cores=2, shift_len=0, save_corrmat=True,)
		
		
	def waveform_correlation(self, input_dir, station_list, lp, hp, t_before, t_after, 
							max_lag, cc_plots=False):



		if isinstance(input_dir, str) == False:
			input_dir = str(input)

		else:
			pass


		for station in station_list:
			# Opening output_file with correlations
			Ofile_name = "xcorrdata_%s.dat" %station

			if os.path.exists(Ofile_name):
				os.system("rm %s" %Ofile_name)
			else:
				print("Creating xcorrdata.dat for station %s" %station)

			Ofile = open(Ofile_name, 'a')

			print('#<----------------------------------------->#')
			print('Working CCs with data from Station = %s' %station)

			sta_str =  "%s*.sac" %station
			data_path = os.path.join(input_dir, station, sta_str)
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

			for jj, uni in enumerate(unique):
				position = tuple(uni)

				event1_, event2_ = uni[0], uni[1]
				selid1, selid2 = Event_id[uni[0]], Event_id[uni[1]]
				event1, event2 = st[event1_], st[event2_]

				wave1, ptime1, rtime1 = get_ptime(event1)
				wave2, ptime2, rtime2 = get_ptime(event2)

				rtime1 = rtime1.strftime('%b-%dT%H:%M,%Y')
				rtime2 = rtime2.strftime('%b-%dT%H:%M,%Y')

				Ofile.write("#   %s   %s  0.000\n" %(selid1, selid2))
				tau, coeff = xcorr_pick_correction(ptime1, wave1, ptime2, wave2, 
														t_before=t_before,
														t_after=t_after, 
														cc_maxlag=max_lag, 
														filter="bandpass",
														filter_options={'freqmin': lp, 'freqmax': hp},
														plot=False)
				CCs.append(coeff)
				Taus.append(tau)
				Times1.append(rtime1)
				Times2.append(rtime2)
				X.append(event1_)
				Y.append(event2_)

				Ofile.write('{:>8} {:>8} {:>8} {:>3}\n'.format(station, np.around(tau, 3), np.around(coeff, 3), "P"))
				Ofile.write("OV\n")

				print("Event pair: %s <---> %s, CC: %s, tau: %s, Times: %s <---> %s" %(selid1, selid2, 
										f'{np.around(coeff, 3):.3f}', f'{np.around(tau, 3):.3f}', 
										rtime1, rtime2,))

			Ofile.close()

		unify_xcorr_tables(station_list=station_list)

	# 	#preparing the event list/catalog in the correct format for running GrowClust
		_prep_in_evtlist(seismic_catalog=self.event_locations, min_mag=self.min_mag, max_mag=self.max_mag)

		
if __name__ == '__main__':
    print("Running Waveform Classifier!")
