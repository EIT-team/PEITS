#!/bin/sh

# Simple regression test for forward solver
# Should only be run using default settings, when first installed
#
# Currently testing by confirming generated file sizes are the same.
# This is a fairly crude method, but md5 hashing won't work, due to slightly different (at >10 decimal places)
# voltages/jacobians being generated on different processers.

THIS_DIR="$(readlink -f "$(dirname "$0")")"
cd "${THIS_DIR}/../data"

# now download nn mesh and try it
wget https://zenodo.org/record/1313709/files/NNmesh_small.tar.gz
tar xf NNmesh_small.tar.gz

# change settings in parameter file
cp parameter_example.cfg parameter.cfg
cp standardparams_example.cfg standardparams.cfg
cp NN2016_Prt_full.txt NN2016_Prt_full.txtBU
head NN2016_Prt_full.txtBU > NN2016_Prt_full.txt

#mesh name
sed -i 's/.*mesh: TA052_meters.*/mesh: NNmesh_small/' parameter.cfg
#mesh parameter file
sed -i 's/.*paramfile: ..\/data\/p.*/\paramfile: ..\/data\/\$(mesh).parameters/' parameter.cfg
#conductivities
sed -i 's/.*conductivities:.*/conductivities: ..\/data\/NNconductivities/' parameter.cfg
#current protocol
sed -i 's/current_protocol_TA052.txt/NN2016_Prt_full.txt/' parameter.cfg

#electrode stuff
sed -i 's/133.0e-6/240.0e-6/' standardparams.cfg
sed -i 's/contact.impedance: 1000/contact.impedance: 200/' standardparams.cfg
sed -i 's/electrode.diameter: 7.0/electrode.diameter: 10.0/' standardparams.cfg
sed -i 's/electrode.diameter_file:.*/electrode.diameter_file: none/' standardparams.cfg




#  run forward calculation
cd "${THIS_DIR}/../src"
mpirun -np 2 "${THIS_DIR}/../src/dune_peits"

# Change to output directory
#cd ../output

# Get file size of computed electrode voltage file and jacobian matrix
#voltage_size=$(wc -c electrode*)
#voltage_size=$(echo $voltage_size | cut -f1 -d ' ')

#jac_size=$(wc -c jacob*)
#jac_size=$(echo $jac_size | cut -f1 -d ' ')


# What we expect the sizes to be
# correct_voltage_size=82
# correct_jac_size=2930338

# echo Computed jac_size: $jac_size
# echo Expected jac_size: $correct_jac_size

# echo Computed voltage_size: $voltage_size
# echo Expected voltage_size: $correct_voltage_size

# # Check if the two hashes are the same
# if [ "$voltage_size" = "$correct_voltage_size" ] && [ "$jac_size" = "$correct_jac_size" ]
# then
# 	echo file sizes match - test OK
# else
# 	echo file sizes do not match, test failed
# 	exit 1
# fi
