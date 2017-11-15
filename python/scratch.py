import os
import glob
import subprocess
import random
from run_many_PEITS import *


def main():
	'''direc = ('../data/01001/')
	mesh_names = get_mesh_names(direc)
	
	pert_types = [{"pathology" : 'haemo', "conductivity" : 1}, {"pathology" : 'isch', "conductivity" : -1}]
	pert_sizes = [5e-3, 10e-3, 15e-3]

	for mesh in mesh_names:
		mesh_path = direc + mesh

		pert_info = {}
		pert_info["do_perturbation"] = random.randint(0, 1)
	
		if pert_info["do_perturbation"]:
			pert_info["location"] = get_pert_location(mesh_path)
			pert_info["size"] = random.choice(pert_sizes)
			pert_info["type"] = random.choice(pert_types)
		print pert_info'''

	data_dir = './data/'

	min_spec_error = -0.1
	max_spec_error = 0.1
	tissue_spectral_error = np.random.uniform(min_spec_error, max_spec_error, size=(7,1))
	conductivity = generate_conductivities	(tissue_spectral_error)
	
	print tissue_spectral_error
	for c in conductivity:
		print c

	for freq in xrange(12):
		 
		write_conductivity_spectra(data_dir, freq, conductivity)



if __name__ == "__main__":
	main()
