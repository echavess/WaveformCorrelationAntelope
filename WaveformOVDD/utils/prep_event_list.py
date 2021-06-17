''' 
This code is part of the Waveform_relocation code for Double difference earthquake relocation
Author: Dr. Esteban J. Chaves
Volcanological and Seismological Observatory of Costa Rica
OVSICORI-UNA
2019
'''


import numpy as np 
import datetime
import os, shutil


def _prep_in_evtlist(seismic_catalog, min_mag, max_mag):

	Day, Time, Latitude, Longitude, Depth, Magnitude, evid = np.genfromtxt(seismic_catalog, unpack=True,
																		 usecols=(0, 1, 2, 3, 4, 5, 9,),
																		 dtype='str')

	Magnitude = Magnitude.astype(np.float)
	Latitude = Latitude.astype(np.float)
	Longitude = Longitude.astype(np.float)
	Depth = Depth.astype(np.float)

	Time_vec = []

	for i, day in enumerate(Day):
		temp = day + "T" + Time[i]
		Time_vec.append(datetime.datetime.strptime(temp, "%Y%m%dT%H%M%S%f"))


	Nmag = []
	Nlat =[]
	Nlon =[]
	Ndep = []
	nTime = []
	nevid = []

	# Selecting events by magnitudes for extracting templates
	for i in range(len(Magnitude)):

		if Magnitude[i] >= min_mag and Magnitude[i] <= max_mag:
			Nmag.append(Magnitude[i])
			Nlat.append(Latitude[i])
			Nlon.append(Longitude[i])
			Ndep.append(Depth[i])
			nTime.append(Time_vec[i])
			nevid.append(evid[i])

	text_file_out = 'evlist.dat'
	out_path = os.path.join('IN', text_file_out)

	if os.path.exists(text_file_out):
		os.remove(text_file_out)

	if os.path.exists(out_path):
		os.remove(out_path)


	outfile = open(text_file_out,'a')

	for gg, t in enumerate(nTime):

		msec = (t.microsecond//1000)

		_msec_str = str(msec)
		l = len(_msec_str)

		if l < 3:
			if l ==1:
				msec = msec * int(100)
			elif l  == 2:
				msec = msec * int(10)
		else:
			msec = msec

						
		outfile.write("{:4d}  {:2d} {:2d}  {:2d} {:2d}  {:2d}.{:02d}  {:06.4f} {:06.4f}  {:.3f}    {:.2f}  0.000  0.000  0.000    {:10.8}\n".format(t.year, t.month, t.day, t.hour, t.minute, t.second, 
																					  msec, Nlat[gg], Nlon[gg], Ndep[gg], Nmag[gg], nevid[gg]))

	outfile.close()
	#shutil.move(text_file_out, out_path)

 
