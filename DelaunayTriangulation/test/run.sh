make
time ./test < $1 > $2
python3 ../vis/vis.py $1 $2 $3
make clean
