# Unificar las tablas de correlación para crear el archivo xcorrdata.dat
# El nombre de los archivos de entrada es el siguiente:
# xcorrdata_%station.dat
import os

def unify_xcorr_tables(station_list):

	def count_events(file):
		count = 0
		with open(_file_) as f:
			for line in f:
				if line.startswith("#"):
					count +=1

			return count

	def all_same(items):
		return all(x == items[0] for x in items)

	N = []
	for station in station_list:
		_file_ = "xcorrdata_%s.dat" %station
		N.append(count_events(_file_))

	station_N = station_list[0]

	if all_same(N) is True:
		print("same N events for all stations")
		station_N = station_list[0]

	else:
		max_n_idx = N.index(max(N))
		station_N = station_list[max_n_idx]
		print("station %s contains more events than the rest" %station_N)

	# Reading main file
	_main_file = "xcorrdata_%s.dat" %station_N

	ev1 = []
	ev2 = []

	with open(_main_file) as f:
		for line in f:
			if line.startswith("#"):
				events = line.split(" ")
				ev1.append(events[3])
				ev2.append(events[6])

	# Creating a new -final output file
	out_file_name_main = 'xcorrdata.dat'
	ofinal = open(out_file_name_main, 'a')
	print("Preparing xcorrdata file!")

	for ii, e in enumerate(ev1):
		ofinal.write("#   %s   %s  0.000\n" %(e, ev2[ii]))

		for ksta in station_list:
			file_corr_results = "xcorrdata_%s.dat"%ksta

			if not os.path.exists:
				_msg2 = "There is no file with correlation results for station %s " %ksta
				raise IOError(_msg2)

			with open(file_corr_results) as f:
				for line in f:
					
					_line_a = "#   %s   %s  0.000" %(e, ev2[ii])
					_line_b = "#   %s   %s  0.000" %(ev2[ii], e)
					
					if line.startswith(_line_a) or line.startswith(_line_b):	
						for line in f:
							events = line.split(' ')
							events = list(filter(None, events))

							if len(events) >= 4:
								ofinal.write("    {:5.5}  {:7.6}  {:8.7} P\n".format(events[0], events[1], events[2]))

							if line.startswith("OV"):
								break

	ofinal.close()

