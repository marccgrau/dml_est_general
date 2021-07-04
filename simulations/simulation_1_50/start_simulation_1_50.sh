#!/bin/bash

echo part_1_sim1_50.R
chmod u+rwx part_1_sim1_50.R

echo part_ff_sim1_50.R
chmod u+rwx part_ff_sim1_50.R

Rscript ./part_1_sim1_50.R
echo "Congratulations! Part 1 of simulation 1 executed successfully. Starting with part 2."

for i in {2..11..1}
	do
    	Rscript ./part_ff_sim1_50.R
	
	if [[ $i -lt 10 ]]
	then
	echo "Congratulations! Part $i of simulation 1 executed successfully."
	else 
	echo "Congratulations! Part $i of simulation 1 executed successfully. The simulation ended successfully."
	fi
done