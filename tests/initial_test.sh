# Simple regression test for forward solver
# Should only be run using default settings, when first installed

# Empty output directory
cd ../output
rm *.*

#  run forward calculation
cd ../src
mpirun -np 2 ./dune_peits

cd ../output

# Calcuate hash of computed electrode voltage file
md5hash=$(md5sum electrode*)
md5hash=$(echo $md5hash | cut -f1 -d ' ')

# What we expect the hash to be
correct_hash=acdcbf9cc769dca75fc983aba6f6067a

echo Computed hash: $md5hash
echo Expected hash: $correct_hash

# Check if the two hashes are the same
if [ "$md5hash" = "$correct_hash" ]
then
	echo Hashes match - test OK
else
	echo Hashes do not match, test failed
fi
