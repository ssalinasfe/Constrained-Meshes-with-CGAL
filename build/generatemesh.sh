point_file="LR$1.node"
node_file="LR$1.1.node"
ele_file="LR$1.1.ele"
neigh_file="LR$1.1.neigh"
output_off_triangle="LR$1.1.off"
output="LR$1.2"
output_off_polylla="LR$1.2.off"

echo -n "Generating points...\n"
cd ../data
python3 LRandom.py $1
cd ../build


echo -n "Generating triangulation...\n"
make && ./triangulation ../data/${point_file} "../data/LR$1"


echo -n "Generating mesh...\n"
./Polylla ../data/${node_file} ../data/${ele_file} ../data/${neigh_file} ../data/${output} && geomview -nopanels ../data/${output_off_polylla}
echo "done"
