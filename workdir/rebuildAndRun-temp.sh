cd ..
build
cd workdir
analysis -a Toy -f data/Ntuple.root -s test -c collections.txt -n -1000 >& Toy7.log &
