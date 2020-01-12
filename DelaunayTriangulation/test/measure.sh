make measure

for x in {100,200,400,800,1600,3200,6400,12800,25600,51200,102400};
do
	for i in $(seq 1 1 100);
	do
		echo $x $i
		python3 ../gen/gen.py $x 0 1000 > test_case
		./measure < test_case >> log
	done

	cat log | awk '{m+=$1} END{print m/NR;}' >> time-result.txt

	rm log
done

rm test_case
make clean
