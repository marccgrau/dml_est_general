#!/bin/bash

echo part_1_sim2_10.R
chmod a+rwx part_1_sim2_10.R

echo part_ff_sim2_10.R
chmod a+rwx part_ff_sim2_10.R

Rscript ./part_1_sim2_10.R
echo "Congratulations! Part 1 of simulation 2 executed successfully. Starting with part 2."

for i in {2..10..1}
	do
    	Rscript ./part_ff_sim2_10.R
	
	if [[ "$i" == "10" ]]
	then
	echo "Congratulations! Part $i of simulation 2 executed successfully. The simulation ended successfully."
	else 
	echo "Congratulations! Part $i of simulation 2 executed successfully."
	fi
done