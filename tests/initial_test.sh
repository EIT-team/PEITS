cd ../output
rm *.*

cd ../src
mpirun -np 2 ./dune_peits

cd ../output

md5hash=$(md5sum electrode*)
md5hash=$(echo $md5hash | cut -f1 -d ' ')

correct_hash=acdcbf9cc769dca75fc983aba6f6067a

echo Computed hash: $md5hash
echo Expected hash: $correct_hash


if [ "$md5hash" = "$correct_hash" ]
then
	echo Hashes match - test OK
else
	echo Hashes do not match, test failed
fi
