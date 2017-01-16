# Simple regression test for forward solver
# Should only be run using default settings, when first installed
#
# Currently testing by confirming generated file sizes are the same.
# This is a fairly crude method, but md5 hashing won't work, due to slightly different (at >10 decimal places)
# voltages/jacobians being generated on different processers.

#  run forward calculation
cd ../src
mpirun -np 2 ./dune_peits

# Change to output directory	
cd ../output

# Get file size of computed electrode voltage file and jacobian matrix
voltage_size=$(wc -c electrode*)
voltage_size=$(echo $voltage_size | cut -f1 -d ' ')

jac_size=$(wc -c jacob*)
jac_size=$(echo $jac_size | cut -f1 -d ' ')


# What we expect the sizes to be
correct_voltage_size=82
correct_jac_size=2930338

echo Computed jac_size: $jac_size
echo Expected jac_size: $correct_jac_size

echo Computed voltage_size: $voltage_size
echo Expected voltage_size: $correct_voltage_size

# Check if the two hashes are the same
if [ "$voltage_size" = "$correct_voltage_size" ] && [ "$jac_size" = "$correct_jac_size" ]
then
	echo file sizes match - test OK
else
	echo file sizes do not match, test failed
fi
