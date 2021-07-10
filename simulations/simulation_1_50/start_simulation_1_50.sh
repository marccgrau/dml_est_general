#!/bin/bash

echo part_1_sim1_50.R
chmod a+rwx part_1_sim1_50.R

echo part_ff_sim1_50.R
chmod a+rwx part_ff_sim1_50.R

Rscript ./part_1_sim1_50.R
echo "Congratulations! Part 1 of simulation 1 executed successfully. Starting with part 2."

for i in {2..10..1}
	do
    	Rscript ./part_ff_sim1_50.R
	
	if [[ "$i" == "10" ]]
	then
	echo "Congratulations! Part $i of simulation 1 executed successfully. The simulation ended successfully."
	else 
	echo "Congratulations! Part $i of simulation 1 executed successfully."
	fi
done