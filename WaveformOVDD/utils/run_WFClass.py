import numpy as np
import os
from WaveformOVDD.utils.WaveformOVDD import Waveform_relocation
from WaveformOVDD.utils.prep_event_list import _prep_in_evtlist


stations = ['LAFE', 'JACO', 'ORTG']
Catalogs = ['entrada_golfo.event.cat', 'entrada_golfo.absolute.data']
min_mag=4.0
max_mag=5.0

step1=False
step2=True


WFC = Waveform_relocation(Miniseed_dir='miniseed', Seismic_Catalogs=Catalogs, min_mag=min_mag, max_mag=max_mag)

# ''' 
# 1. Extracting the waveforms from the miniseed 24-hr files using the catalogs of events
# output1: Sac data organized in magnitude bins for each station in list. 
# output2: list with event ids for each station
# ''' 

if step1 is True:

	for ista in stations:

		WFC.get_templates_from_catalog(station=ista, stime=5, endtime=70, save_waves=False)

else:
	print("stations were downloaded already...Checking the P-wave picks")

'''
2. If SAC data is not picked yet, please comment the next function. 
Make sure your data is picked (P-Wave) before running next. 
'''
# 2. Correlating the vertical component of the seismograms for all events in catalog and generating the input
# files for GrowClust

if step2 is True:
	directory = "Waveforms_%s_%s" %(min_mag, max_mag)
	WFC.waveform_correlation(input_dir=directory, station_list=stations, 
							lp=1, hp=10.0, 
							t_before=0.2, 
							t_after=0.6, 
							max_lag=1.0)

# Running GrowClust
#_cmd_ = 'growclust growclust.inp'
#os.system(_cmd_)


