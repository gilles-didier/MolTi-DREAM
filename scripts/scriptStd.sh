#usage ./script.sh <partition in clust n see format> <graphs directory>
#write the new partition in <partition in clust n see format>_new.cls, still in ClustNsee format
#example : ./script.sh outFile_p5.cls ../datasets/BiologicalNetworksMay2015
#you have to set the path to the executable according to your working directory

# -m set the size of the classes we have to split
../src/extract_part -m 100 -f c $1 `basename $1 | cut -d'.' -f1`
for f in *.tmp; do
	echo $f
	for g in $2/*; do
		h=`basename $g | cut -d'.' -f1`'_'`echo $f | cut -d'.' -f1`.new
		../src/extract_graph $g $f $h
	done
	h=`ls *.new`
# here one can change the molti-console option (resolution parameter etc.)
	../src/molti-console -r 0 -p 1 -f s -o `basename $f | cut -d'.' -f1`.clstmp $h
	rm *.new
	rm *.csv
done
rm *.tmp
cat *.clstmp > `basename $1 | cut -d'.' -f1`'_'tmp.cls
# 
rm *.clstmp
# back to the ClustNsee format
../src/stdTocNs `basename $1 | cut -d'.' -f1`'_'tmp.cls `basename $1 | cut -d'.' -f1`'_'new.cls
rm `basename $1 | cut -d'.' -f1`'_'tmp.cls
