dirs=(10v 3t 5s 9t 2m 4r 6s)

for D in ${dirs[*]}
do
	cd $D
	qsub submit_gtxtitan
	cd ../
done
