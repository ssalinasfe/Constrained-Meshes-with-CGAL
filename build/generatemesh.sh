point_file="square_$1.node"
node_file="square_$1.1.node"
ele_file="square_$1.1.ele"
neigh_file="square_$1.1.neigh"
output_off_triangle="square_$1.1.off"
output="square_$1.2"
output_off_polylla="square_$1.2.off"

echo -n "Generating points...\n"
cd ../data
#python3 LRandom.py $1
python3 10000x10000RandomPoints.py $1
cd ../build


#echo -n "Generating triangulation...\n"
make && ./triangulation ../data/${point_file} "../data/square_$1"


#echo -n "Generating mesh...\n"
./Polylla ../data/${node_file} ../data/${ele_file} ../data/${neigh_file} ../data/${output} 
#echo "done"

#geomview -nopanels ../data/${output_off_polylla}
