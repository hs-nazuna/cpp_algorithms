make verify

for i in $(seq 1 1 100):
do
	python3 ../gen/gen.py 1000 0 100 > test_data.in
	./verify < test_data.in
done

rm test_data.in

make clean
