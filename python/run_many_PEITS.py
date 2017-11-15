import os
import glob
import subprocess
import random
import json
import numpy as np
import matplotlib.pyplot as plt
import optparse

def get_mesh_names(directory):

	dir_contents =  os.listdir(directory)
	mesh_names = []

	for file in dir_contents:
		if file.endswith(".dgf"):
			mesh_names.append(file[:-4])

	return mesh_names


def update_parameter_file(mesh_path, pert_info, freq):

	param_file = './data/param_deform'
	f = open(param_file, 'w')


	f.write('fem.io.macroGrid: .'+ mesh_path +'.dgf\n')
	f.write('electrode.use_node_assignment: false\n')
	f.write('electrode.positions: .' + mesh_path + '.electrodes\n')
	f.write('fem.assign_conductivities: true\n')

	# First time through need to create paritions
	if pert_info["num"] == 0:
		f.write('fem.io.loadPartitions: false\n')
	# The rest of the time we can load them from the first run
	else:
		f.write('fem.io.loadPartitions: true\n')

	if pert_info["do_perturbation"]:
		f.write('mesh.perturbation: true\n')
	else:
		f.write('mesh.perturbation: false\n')

	# true for multiplication with value, false for setting absolte value
	f.write('mesh.perturbation.multORabs: false\n')
	f.write('mesh.perturbation.value: ' + str( pert_info["type"]["conductivity"][freq]) + '\n')
	f.write('mesh.perturbation.radius: ' + str( pert_info["size"] ) + '\n' )

	x,y,z = pert_info["location"];

	f.write('mesh.perturbation.pos_x: ' + str(x) + '\n')
	f.write('mesh.perturbation.pos_y: ' + str(y) + '\n')
	f.write('mesh.perturbation.pos_z: ' + str(z) + '\n')



	#Append existing parameters to the file
	existing_parameter_file = open(mesh_path + '.parameters', 'r')
	f.write(existing_parameter_file.read())
	existing_parameter_file.close()

	f.close()

	#f = open(param_file, 'r')
	#print f.read()

def get_newest_file(path, glob_pattern):
#Find the newest file in a directory that matches glob_pattern
	#print "Newest file: "
	newest = min(glob.glob(path + glob_pattern), key = os.path.getctime)
	return newest


def rename_newest_file(path, new_name):
	#Find the newest electrode voltage file and rename it, to contain the mesh name
	#and delete the sigma vector file
	print "Newest files: "
	newest = min(glob.glob(path + 'elec*.bin'), key = os.path.getctime)
	print newest
	print path + new_name
	os.rename(newest, path + new_name)

	#Delete sigma files
	for f in glob.glob(path + 'sigma*'):
		os.remove(f)


def get_pert_location(mesh_name):
	cells_fname = mesh_name + '.cell_centres';
	n_cells = file_len(cells_fname)

	found_brain = False

	#Data is stored as string
	brain_ind = ['1','2']

	rand_cell = random.randint(1, n_cells)
	line = read_line(cells_fname, rand_cell)
	cell_data =  line.split()
	while not found_brain:
		rand_cell = random.randint(1, n_cells)
		line = read_line(cells_fname, rand_cell)
		# Line of cell_centres should be x y z type
		# We want the type to be brain
		cell_data =  line.split()

		if cell_data[-1] in brain_ind:
			found_brain = True

	x = float(cell_data[0])
	y = float(cell_data[1])
	z = float(cell_data[2])

	x,y,z = convert_to_metres(x,y,z)

	return [x,y,z]


def convert_to_metres(x,y,z):
# Units of distance will either be in mm or metres
# if in mm, convert to metres
	if (x > 1):
		print "Seems units are in mm, converting to metres"
		return x/1000, y/1000, z/1000

	return x,y,z


def get_pert_spectrum(spectral_error):
		# Add a spectral error to the perturbation conductivity
		blood = [0.6972] * 12
		ischaemia = [ 0.022, 0.0277, 0.0509, 0.0758, 0.0837, 0.0897, 0.0939, 0.0972, 0.0983, 0.0983, 0.1049, 0.1096]

		blood = [x + spectral_error * x for x in blood]
		ischaemia = [x + spectral_error * x for x in ischaemia]

		return blood, ischaemia


def get_pert_info(mesh_path, spectral_error):

	blood, ischaemia = get_pert_spectrum(spectral_error)

	pert_types = [	{"pathology" : 'haemo', "conductivity" : blood},
			{"pathology" : 'isch', "conductivity" : ischaemia}]

	pert_sizes = [10e-3, 15e-3]

	pert_info = {}
	#TODO: choose when to do a perturbation
	pert_info["do_perturbation"] = 1#random.randint(0, 1)

	pert_info["location"] = get_pert_location(mesh_path)
	pert_info["size"] = random.choice(pert_sizes)
	pert_info["type"] = random.choice(pert_types)

	return pert_info


def read_line(fname, line_num):

	# Open file and go through lines until we reach the one we want
 	with open(fname) as f:
		for i, l in enumerate(f):
			if i == line_num:
				return l

	print "Line not found in file"
	return -1


def file_len(fname):
	#print fname
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1


def generate_conductivities(spectral_error):
	#Returns a set of condutivity values for different issues, adjusted for spectral error

	# These are the base values (not sure where the origin was, but we normally use these values)
	cond = []
	cond.append( [0.0921, 0.0945, 0.1009, 0.1121, 0.1140, 0.1140, 0.1155, 0.1212, 0.1234, 0.1234, 0.1311, 0.1437] )# White
	cond.append ([2 * x for x in cond[0]]) # Grey = double white matter
	cond.append ([ 1.79] * 12) # CSF
	cond.append ([3 * x for x in cond[0]] )# Dura =  triple white matter
	cond.append ([0.018] * 12 )# Skull
	cond.append ([0.0001] * 12) # Air
	cond.append ([0.225, 0.225, 0.23, 0.25, 0.265, 0.28, 0.35, 0.39, 0.405, 0.41, 0.42, 0.425] )# Scalp

	#Adjust for spectral error
	cond = cond + cond * spectral_error

	return cond

def write_conductivity_spectra(data_dir, freq, cond):
	# Write the conducitvity data, for a particular frequency to the file needed by PEITS
	# data_dir : where to write the file
	# cond: conductivity values for all tissue
	# freq: the frequency of interest, used as index to cond array
	cond_file = data_dir + 'conductivities'

	with open(cond_file, 'w') as f:
		for i in cond:
			f.write(str(i[freq]) + '\n')


#https://stackoverflow.com/questions/12517451/automatically-creating-directories-with-file-output
def mkdir_p(path):
# Check if a directory exists. If it doesn't create it.
	if not os.path.exists(path):

	    try:
        	os.makedirs(path)
	    except OSError as exc: # Python >2.5
        	if exc.errno == errno.EEXIST and os.path.isdir(path):
 			raise


def load_electrode_voltages_binary(fname):
	#Reimplemented from MATLAB function with same name
	# Read in data from PEITS output .bin file.

        #Check that magicstr, int and double written correctly by forward solver
        # then load data
        with open(fname, 'rb') as f:
            c = np.fromfile(f, np.uint8, count=3)
            magicstr = ''.join([str(unichr(x)) for x in c])

            if not magicstr == 'DEV':
                print ('Read magicstr doe not indicate Dune Sparse Matrix!')

            magicint = np.fromfile(f,np.int32, count=1)
            magicdouble = np.fromfile(f,np.float,count=1)

            if not magicint == 111 or not magicdouble == 111.0:
                print "Magic nubmers not read correctly, check the binary format in this reading routing!!"

            ncols = np.fromfile(f,np.int32, count=1)
            nrows = np.fromfile(f,np.int32, count=1)

            v = np.fromfile(f, np.float, count = (nrows*ncols))
            v.reshape(nrows,ncols)

            c = np.fromfile(f, np.uint8, count=3)
            eofstr = ''.join([str(unichr(x)) for x in c])

            if not eofstr == 'EOF':
                print "Read eofstr does not indicate end of binary file!"

        return v

def save_spectrum_image(output_name, tissue_spectrum, pert):
	#Plot data and save to file
	plt.semilogy(tissue_spectrum.T)
	plt.semilogy(pert["type"]["conductivity"], linestyle = '--')

	plt.legend([	"White Matter", "Grey Matter", "CSF",
	 			"Dura", "Skull", "Air", "Scalp",
				 pert["type"]["pathology"]],
				 loc = 'lower right')

	plt.ylim([1e-5, 1e1])
	plt.savefig(output_name + '.jpg')
	plt.close()

def check_valid_voltages(voltages):
	max_voltage = 0.5 # Maximum sensible voltage value
	if max(abs(voltages)) > max_voltage:
		return 0
	else:
		return 1

def main():

	#Parse commandline arguments
	parser = optparse.OptionParser()
	parser.add_option(	'--verbose', dest="verbose", action="store_true", default=False,
						help="Don't print PEITS output")

	(options, args) = parser.parse_args()

	#Load list of all the meshes
	#We want to run multiple simulations on each mesh

	output_dir = './output/'
	data_dir = './data/'
	direc = ('./data/meshes/')
	mesh_names = get_mesh_names(direc)

	num_perturbations = 1
	num_freqs = 12

	sim_failures = [] #Info on simulations that produce 'bad' output (fails check_valid_voltages())
	successful_sims = 0

	# Simulate perturbations in each mesh, at 12 frequencies
	for mesh in mesh_names:
		print("Mesh: " + str(mesh))
		mesh_path = direc + mesh

		for n in xrange(num_perturbations):

			print ("Perturbation number: %s" % n)

			# Output file name to save voltages and sim info
	 		output_name = output_dir + mesh + "_pert_" + str(n)

			voltages = []

			# Generate a spectral error between -10% and +10% for each of the 7 'normal' tisuses
			# Need to use size=(7,1) rather than just size=7 to allow numpy to broadcast the array
			min_spec_error = -0.1
			max_spec_error = 0.1
			tissue_spectral_error = np.random.uniform(min_spec_error, max_spec_error, size=(7,1))
			conductivity = generate_conductivities(tissue_spectral_error)

			# Pick perturbation size, location and type
			pert_spectral_error = np.random.uniform(min_spec_error, max_spec_error)
			pert_info = get_pert_info(mesh_path, pert_spectral_error)
			pert_info["num"] = n
			pert_info["mesh_name"] = mesh

			# If any of the simulated voltages fail validation, set this to false
			# So that we don't save unnecssary values
			voltages_ok = True

			for f in xrange(num_freqs):


				print("Frequency: %s " % f)

				mesh_output_file = mesh_path + "freq_" + str(f)

				update_parameter_file(mesh_path, pert_info, f)
				write_conductivity_spectra(data_dir, f, conductivity)

				#os.system('mpirun -np 12 ./src/dune_peits')
				#TODO: Fix relative path/script call directory


				if options.verbose:
					#Use os module,which prints all output to screen
					os.system('sh ./matlab_c_call.sh')
				else:
					# Use subprocess, which doesn't print anything
					subprocess.check_output(['sh', './matlab_c_call.sh'])

				#Read voltages from file
				voltage_bin_file = get_newest_file(output_dir, 'elec*.bin')


				#Delete sigmavector files, don't need them
				for filename in glob.glob(output_dir + 'sigma*'):
					os.remove(filename)

				# Check if the voltages are sensisble,
				# if not don't continue with this set of simulations, exit for loop
				# If they are sensible, append them to the voltages list
				v = load_electrode_voltages_binary( voltage_bin_file)
				#Delete binary voltage file
				os.remove(voltage_bin_file)

				if not check_valid_voltages(v):
					print "Simulated voltages did not pass sanity check, discarding"
					voltages_ok = False
					# Store some data about this simulation
					sim_failures.append([mesh, f, pert_info])

					break

				voltages.append(v)

			#Only save the data if all of the simulated voltages are good
			if voltages_ok:

				#Save voltages as numpy array and write info about the perturbation for later reference
				voltages_numpy = np.array(voltages)
				print "Saving voltages as numpy array."
				numpy_output_file = output_name + '.npy'
				sim_info_file = output_name + '.txt'

				np.save(numpy_output_file, voltages_numpy) # Save voltages
				with open(sim_info_file, 'w') as file: # Save sim info
					file.write(json.dumps(pert_info))
				# Save conducivity spectrum plot as a jpg
				save_spectrum_image(output_name, conductivity, pert_info)
				successful_sims = successful_sims + 1

print("OK simulations: %s	Failed simulations: %s" % (successful_sims, len(sim_failures)))

# Save info on failures to JSON
sim_failures_file = "failures.json"
with open(sim_failure_file, 'w') as file:
	file.write(json.dumps(sim_failures))

	



if __name__ == "__main__":
	main()
